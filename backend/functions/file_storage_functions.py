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
                    "extension": "folder",
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
            if delete_path != Path("/"):
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
    Handles file uploads to Cloud Storage and creates corresponding Firestore documents.
    
    Triggered when a file is uploaded to the 'twocube-web.appspot.com' bucket.
    Handles two types of paths:
        - "uploads/{userID}/{uploadPath}" for user files
        - "sessionImages/{userID}/{imageName}" for session images
    
    For user files, the Firestore document will be created under 'users/{userID}/files/{documentID}'.
    For session images, the Firestore document will be created in the 'sessionImage' collection.
    """
    # Retrieve cloud storage path
    cloud_storage_path = event.data.name

    # Define regex patterns for both paths
    user_file_pattern = r"^uploads/([^/]+)(.+)$"
    session_image_pattern = r"^sessionImages/([^/]+)(.+)$"

    # Match the path with appropriate pattern
    user_file_match = re.match(user_file_pattern, cloud_storage_path)
    session_image_match = re.match(session_image_pattern, cloud_storage_path)

    # Initialize Firestore client
    db = firestore.client()

    if user_file_match:
        # Handle user file upload
        user_uid, upload_path = user_file_match.groups()
        upload_path = Path(upload_path)
        upload_name = upload_path.name
        upload_parent = upload_path.parent.as_posix()

        # Verify user ID exists
        transaction = db.transaction()
        user_doc = db.collection("users").document(user_uid)
        if not user_doc.get().exists:
            raise KeyError(f"User ID {user_uid} does not exist")

        # Check if uploaded file already exists
        file_ref = db.collection("users", user_uid, "files")
        query = (
            file_ref.where(filter=FieldFilter("path", "==", upload_parent))
            .where(filter=FieldFilter("name", "==", upload_name))
            .limit(1)
        ).stream()

        # If file exists, use the existing DocumentReference; otherwise, create a new one
        for doc in query:
            file_doc = doc.reference
            break
        else:
            file_doc = file_ref.document()

        # Create or update the Firestore document
        file_doc.set(
            {
                "id": file_doc.id,
                "name": upload_name,
                "path": upload_parent,
                "size": int(event.data.size),
                "extension": upload_path.suffix,
                "isFolder": False,
                "storageLink": event.data.name,
                "uploadedOn": firestore.SERVER_TIMESTAMP,
            }
        )

        # Ensure folders exist in Firestore
        ensure_folder_exists(transaction, db, user_uid, upload_parent)

    elif session_image_match:
        # Handle session image upload
        user_uid, upload_name = session_image_match.groups()
        upload_path = Path(f"sessionImages/{user_uid}/{upload_name}")
        upload_parent = upload_path.parent.as_posix()
        upload_name = upload_name.lstrip('/')

        # Verify user ID exists
        user_doc = db.collection("users").document(user_uid)
        if not user_doc.get().exists:
            raise KeyError(f"User ID {user_uid} does not exist")

        # Create a new Firestore document in the 'sessionImage' collection
        image_doc_ref = db.collection("sessionImage").document()
        image_doc_ref.set(
            {
                "id": image_doc_ref.id,
                "name": upload_name,
                "path": upload_parent,
                "size": int(event.data.size),
                "storageLink": cloud_storage_path,
                "uploadedOn": firestore.SERVER_TIMESTAMP,
                "owner": user_uid,
            }
        )

    else:
        print("Path does not match any expected pattern. Exiting function.")