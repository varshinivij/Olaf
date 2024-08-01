from firebase_admin import auth, firestore

from firebase_functions.https_fn import (
    CallableRequest,
    FunctionsErrorCode,
    HttpsError,
    Request,
    on_call,
)
from firebase_functions.storage_fn import (
    CloudEvent,
    StorageObjectData,
    on_object_finalized,
)

from google.cloud import storage
from google.cloud.exceptions import NotFound
from google.cloud.firestore import Client, transactional
from google.cloud.firestore_v1.base_query import FieldFilter
from google.cloud.firestore_v1 import Transaction

from collections import deque
from pathlib import Path
import re


def validate_firebase_id_token(req: Request) -> dict:
    """
    Validates the Firebase ID token passed in a request's Authorization header.
    Returns the decoded token as a dict containing user data. Example returned
    payload:
    {
        "name": str
        "picture": str (link)
        "iss": str (link),
        "aud": str,
        "auth_time": int,
        "user_id": str,
        "sub": str,
        "iat": int,
        "exp": int,
        "email": str,
        "email_verified": bool,
        "firebase": {
            "identities": {
                "google.com": list[str],
                "email": list[str],
            },
            "sign_in_provider": str,
        },
        "uid": str,
    }

    To send the ID token in Angular:
    const headers = new HttpHeaders({
      Authorization: `Bearer ${await this.auth.currentUser?.getIdToken()}`,
    });
    this.http.post(url, { param1: 'folder1' }, { headers });

    Not necessary for HTTP Callables, only pure HTTP functions. Not currently
    used but here for your convenience.
    """
    if not req.headers.get("Authorization", "").startswith("Bearer "):
        raise PermissionError(
            "No Firebase ID token passed as Bearer token in Authorization header."
            " Provide the following HTTP header: Authorization: Bearer <Firebase ID Token>"
        )

    id_token = req.headers["Authorization"].split("Bearer ")[1]
    return auth.verify_id_token(id_token, check_revoked=True)


@transactional
def ensure_folder_exists(
    transaction: Transaction, firestore_client: Client, user_uid: str, folder_path: str
) -> None:
    """
    Given a userID and folder path, ensures that all folders along the path
    have corresponding records in Firestore. If not, creates all the folders
    from top-down.

    For example, if folder_path="/folder1/folder2/folder3" but only folder1
    existed, it will create Firestore records for folder2 and folder3.

    This function is wrapped as a transactional function to ensure atomic DB
    operations: either all DB writes succeed, or they all fail, to prevent
    partial writes such as only creating folder3 in case of error.
    """
    folder_queue = deque([Path(folder_path)])

    file_ref = firestore_client.collection("users", user_uid, "files")

    while folder_queue:
        current_folder = folder_queue.popleft()

        if current_folder == Path("/"):  # reached root dir, confirm it exists
            break

        query = (
            file_ref.where(
                filter=FieldFilter("path", "==", current_folder.parent.as_posix())
            )
            .where(filter=FieldFilter("name", "==", current_folder.name))
            .where(filter=FieldFilter("isFolder", "==", True))
            .limit(1)
        ).stream()

        for doc in query:
            break  # if folder exists, we're done
        else:
            # if it doesn't exist, create it
            file_doc = file_ref.document()
            file_doc.set(
                {
                    "id": file_doc.id,
                    "name": current_folder.name,
                    "path": current_folder.parent.as_posix(),
                    "size": 0,
                    "fileExtension": 'folder',
                    "isFolder": True,
                    "storageLink": f"uploads/{user_uid}{current_folder.as_posix()}",
                    "uploadedOn": firestore.SERVER_TIMESTAMP,
                }
            )

            # add parent to be created next
            folder_queue.append(current_folder.parent)


@on_call(region="us-central1")
def request_user_create_folder(req: CallableRequest) -> None:
    """
    HTTP Callable function that creates a folder with a specified user and path.
    Takes two parameters:
        name: name of the folder to create
              format: "folder3" will create a folder in "/folder1/folder2/folder3"
        path: the path to the parent of where to make the folder
              format: "/folder1/folder2", '/' means root dir

    If folder already exists, does nothing.

    This function is an HTTPS Callable function
    (https://firebase.google.com/docs/functions/callable?gen=2nd),
    which is essentially an HTTP endpoint that automatically implements
    Firebase Auth for you and makes the UID easily accessible. You can use
    httpCallable() in Angular to call it instead of using an HTTP request
    and it makes it much easier to implement.

    If you still want to use a normal HTTP endpoint with Firebase Auth,
    use validate_firebase_id_token() defined above.
    """
    try:
        user_uid = req.auth.uid
        folder_name = req.data["name"]
        folder_parent = req.data["path"]
        if not folder_parent.startswith("/"):
            folder_parent = "/" + folder_parent

    except (ValueError, KeyError):
        raise HttpsError(
            code=FunctionsErrorCode.INVALID_ARGUMENT,
            message="User must be logged in and function must be called with arguments 'name': str and 'path': str",
        )

    db = firestore.client()
    transaction = db.transaction()
    folder_path = (Path(folder_parent) / Path(folder_name)).as_posix()

    ensure_folder_exists(transaction, db, user_uid, folder_path)


@on_call(region="us-central1")
def request_user_delete_path(req: CallableRequest) -> None:
    """
    HTTP Callable function that deletes a specified directory.
    Takes one parameters:
        path: the path to the directory to delete.
              format: "/folder1/folder2", '/' means root dir

    If the directory is a folder, deletes all files within the folder.
    Also deletes all corresponding files in Cloud Storage. (Empty folders
    may or not may not exist in Cloud Storage due to its flat structure.)
    """
    try:
        user_uid = req.auth.uid
        delete_path = req.data["path"]
        if not delete_path.startswith("/"):
            delete_path = "/" + delete_path

        delete_path = Path(delete_path)

    except (ValueError, KeyError):
        raise HttpsError(
            code=FunctionsErrorCode.INVALID_ARGUMENT,
            message="User must be logged in and function must be called with arguments 'path': str",
        )

    db = firestore.client()
    transaction = db.transaction()
    file_ref = db.collection("users", user_uid, "files")

    storage_client = storage.Client()
    bucket = storage_client.bucket("twocube-web.appspot.com")

    @transactional
    def execute_as_transaction(transaction):
        # find a single file doc with the same path
        query = (
            file_ref.where(
                filter=FieldFilter("path", "==", delete_path.parent.as_posix())
            )
            .where(filter=FieldFilter("name", "==", delete_path.name))
            .limit(1)
        ).stream()

        # retrieve the first (and only) query result
        for doc in query:
            is_folder: bool = doc.get("isFolder")
            storageLink = doc.get("storageLink")

            # delete firestore record
            doc.reference.delete()

            # delete cloud storage
            # NOTE: transactions are only for firestore so deletions may not sync up.
            try:
                # also, if delete fails (file not found) ignore it
                bucket.blob(storageLink).delete()
            except NotFound:
                pass

            break
        else:  # if no queries exist (didn't break) exit function
            # but if the delete path is root there is no firestore doc,
            # so special case (messy code, but we can structure it better
            # later)
            if delete_path != Path('/'):
                return
            else:
                is_folder = True

        if not is_folder:  # if not folder, deleting once was enough
            return

        # delete firestore + cloud storage for all file docs with a path
        # starting with delete dir

        # see here for how i make this query:
        # https://stackoverflow.com/questions/46573804/firestore-query-documents-startswith-a-string
        lex_next_string = delete_path.as_posix()
        lex_next_string = lex_next_string[:-1] + chr(ord(lex_next_string[-1]) + 1)

        query = (
            file_ref.where(
                filter=FieldFilter("path", ">=", delete_path.as_posix())
            ).where(filter=FieldFilter("path", "<", lex_next_string))
        ).stream()

        for doc in query:
            storageLink = doc.get("storageLink")
            doc.reference.delete()
            try:
                bucket.blob(storageLink).delete()
            except NotFound:
                pass

    execute_as_transaction(transaction)


@on_object_finalized(bucket="twocube-web.appspot.com")
def handle_user_file_upload(event: CloudEvent[StorageObjectData]) -> None:
    """
    Receives a user-uploaded Cloud Storage file and creates a Firestore
    document under users/{userID}/files/{documentID} accordingly.

    Triggered when a file is uploaded to the twocube-web.appspot.com bucket.
    If upload path doesn't match the pattern "uploads/{userID}/{uploadPath}",
    does nothing.

    The newly created Firestore document will contain a random UID. It will
    have the following properties:
        "id": the file's unique UID
        "name": the name of the file,
        "path": the path to the file on the user's dashboard (up to the parent),
                format '/folder1/folder2/file.txt'. '/' path means root dir.
        "isFolder": whether the document is a folder or not
        "size": the size of the document in bytes
        "storageLink": the path to the file in Cloud Storage
        "uploadedOn": the timestamp of when the file was uploaded

    If the Firestore document for the path already exists (overwriting a file
    in Cloud Storage), it will overwrite the corresponding Firestore record
    instead and overwrite it as well. In this way, we are mimicking Cloud
    Storage's mechanics of overwriting duplicate file names.

    One potentially unsolved edge case to consider is what happens if Cloud
    Storage upload succeeds on the frontend but the Firestore record fails to
    create here. For now, it does nothing and re-uploading the same file to
    Cloud Storage should re-run this function and potentially fix it.
    """
    # retrieve cloud storage path
    cloud_storage_path = event.data.name

    # verify that cloud storage path matches the format "uploads/{userID}/{uploadPath}"
    match = re.match(r"^uploads/([^/]+)(.+)$", cloud_storage_path)
    if match is None:
        return

    # retrieve {userID} and {uploadPath} from the above regex
    user_uid, upload_path = match.groups()
    upload_path = Path(upload_path)

    upload_name = upload_path.name
    upload_parent = upload_path.parent.as_posix()

    # verify user ID exists
    db = firestore.client()
    transaction = db.transaction()
    user_doc = db.collection("users").document(user_uid)
    if not user_doc.get().exists:
        raise KeyError(f"User ID {user_uid} does not exist")

    # check if uploaded file already exists
    file_ref = db.collection("users", user_uid, "files")

    query = (
        file_ref.where(filter=FieldFilter("path", "==", upload_parent))
        .where(filter=FieldFilter("name", "==", upload_name))
        .limit(1)
    ).stream()

    # if file exists with the same name and path, use that DocumentReference.
    # else, create a new file.
    for doc in query:
        file_doc = doc.reference
        break
    else:
        file_doc = file_ref.document()

    file_doc.set(
        {
            "id": file_doc.id,
            "name": upload_name,
            "path": upload_parent,
            "size": int(event.data.size),
            "fileExtension": upload_path.suffix,
            "isFolder": False,
            "storageLink": event.data.name,
            "uploadedOn": firestore.SERVER_TIMESTAMP,
        }
    )

    # after creating document, run the batch job to create folders
    ensure_folder_exists(transaction, db, user_uid, upload_parent)
