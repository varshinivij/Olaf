from pathlib import Path
import re

from firebase_functions.https_fn import CallableRequest, on_call
from firebase_functions.storage_fn import (
    CloudEvent,
    StorageObjectData,
    on_object_finalized,
)

from functions.services.file_storage_service import FileStorageService
from functions.services.user_service import UserService
from functions.utils.validation import ValidationError, expect_values_in_request_data


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

    If you want to use a normal HTTP endpoint with Firebase Auth,
    use validate_firebase_id_token() defined in utils/validation.
    """
    if req.auth is None:
        raise ValidationError("User must be logged in to create a folder")

    user_id = req.auth.uid

    args = expect_values_in_request_data(req, ["name", "path"])
    folder_name: str = args["name"]
    folder_parent: str = args["path"]

    if folder_name.startswith("/"):
        folder_name = folder_name[1:]

    folder_path = Path(folder_parent, folder_name)

    file_storage_service = FileStorageService()
    file_storage_service.ensure_folders_exist(user_id, folder_path)


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
    if req.auth is None:
        raise ValidationError("User must be logged in to create a folder")

    user_id = req.auth.uid

    args = expect_values_in_request_data(req, ["path"])
    delete_path: str = args["path"]

    file_storage_service = FileStorageService()
    file_storage_service.delete_path(user_id, delete_path)


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
    # retrieve cloud storage upload path
    cloud_storage_path = event.data.name

    # match upload path with appropriate handler
    user_file_match = re.match(r"^uploads/([^/]+)(.+)$", cloud_storage_path)
    session_image_match = re.match(r"^sessionImages/([^/]+)(.+)$", cloud_storage_path)

    file_storage_service = FileStorageService()
    user_service = UserService()

    # Handle user file upload
    if user_file_match is not None:
        user_id, upload_path = user_file_match.groups()

        if not user_service.get_user_exists(user_id):
            raise ValidationError(f"User ID {user_id} does not exist")

        file_storage_service.create_firestore_upload_record(
            user_id, upload_path, event.data.size, cloud_storage_path
        )

    # Handle session image upload
    elif session_image_match is not None:
        user_id, upload_path = session_image_match.groups()

        if not user_service.get_user_exists(user_id):
            raise ValidationError(f"User ID {user_id} does not exist")

        file_storage_service.create_session_image_record(
            user_id, upload_path, event.data.size, cloud_storage_path
        )

    else:
        print("Path does not match any expected pattern. Exiting function.")
