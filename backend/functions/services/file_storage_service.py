from collections import deque
from pathlib import Path

from firebase_admin import firestore
from google.cloud import storage
from google.cloud.exceptions import NotFound
from google.cloud.firestore_v1.base_query import FieldFilter

from functions.models.session_image import SessionImage
from functions.models.user_file import UserFile


class FileStorageService:
    """
    Service class containing utility methods for file storage endpoints.

    "File storage" specifically refers to user-uploaded files and folders in
    the Files section of the Twocube dashboard, represented as Firestore
    documents and stored in Cloud Storage.
    """

    def __init__(self):
        self.db = firestore.client()

    def normalize_path(self, path: str | Path) -> Path:
        """
        Normalizes a path string, adds leading slash and removes trailing slashes.
        """
        path = Path(path).as_posix()

        if not path.startswith("/"):
            path = "/" + path
        if path.endswith("/"):
            path = path[:-1]

        return Path(path)

    def create_firestore_upload_record(
        self, user_id: str, upload_path: str | Path, upload_size: int, storage_link: str
    ) -> UserFile:
        """
        Creates a Firestore record for a newly uploaded Cloud Storage file
        (for user uploaded files).

        If a Firestore record already exists for the file (same file path),
        overwrites the record.

        upload_path should include the name of the file at the end (the full path).
        """
        upload_path = self.normalize_path(upload_path)
        upload_name = upload_path.name
        upload_parent = upload_path.parent.as_posix()

        files_ref = self.db.collection("users", user_id, "files")

        # check if upload record already exists
        query = (
            files_ref.where(filter=FieldFilter("path", "==", upload_parent))
            .where(filter=FieldFilter("name", "==", upload_name))
            .limit(1)
        ).get()

        # if record already exists, overwrite that record instead
        if len(query) == 1:
            file_doc = query[0].reference
        else:
            file_doc = files_ref.document()

        file_data: UserFile = {
            "id": file_doc.id,
            "name": upload_name,
            "path": upload_parent,
            "size": upload_size,
            "extension": upload_path.suffix,
            "isFolder": False,
            "storageLink": storage_link,
            "uploadedOn": firestore.firestore.SERVER_TIMESTAMP,  # type: ignore
        }

        file_doc.set(file_data)  # type: ignore
        self.ensure_folders_exist(user_id, upload_parent)

        return file_data

    def create_session_image_record(
        self, user_id: str, upload_path: str | Path, upload_size: int, storage_link: str
    ) -> SessionImage:
        """
        Creates a Firestore record for a newly uploaded Cloud Storage file
        (for session image files).

        upload_path should include the name of the file at the end (the full path).
        """
        upload_path = self.normalize_path(upload_path)
        upload_name = upload_path.name
        upload_parent = upload_path.parent.as_posix()

        image_ref = self.db.collection("sessionImage").document()

        image_data: SessionImage = {
            "id": image_ref.id,
            "name": upload_name,
            "path": upload_parent,
            "size": upload_size,
            "storageLink": storage_link,
            "uploadedOn": firestore.firestore.SERVER_TIMESTAMP,  # type: ignore
            "user_id": user_id,
        }

        image_ref.set(image_data)  # type: ignore

        return image_data

    def ensure_folders_exist(self, user_id: str, folder_path: str | Path) -> None:
        """
        Given a userID and folder path, ensures all folders along the path
        have corresponding records in Firestore. If not, creates all folders.

        For example, if folder_path="/folder1/folder2/folder3" but only folder1
        exists, it will create Firestore records for folder2 and folder3.
        """
        files_ref = self.db.collection("users", user_id, "files")

        count = 0
        batch = self.db.batch()

        folder_queue = deque([self.normalize_path(folder_path)])

        while folder_queue:
            current_folder = folder_queue.popleft()
            if current_folder == Path("/"):  # reached root dir, we're done
                break

            query = (
                files_ref.where(
                    filter=FieldFilter("path", "==", current_folder.parent.as_posix())
                )
                .where(filter=FieldFilter("name", "==", current_folder.name))
                .where(filter=FieldFilter("isFolder", "==", True))
                .limit(1)
                .get()
            )

            if len(query) == 1:  # path exists, we're done
                break

            # path doesn't exist, create folder
            file_doc = files_ref.document()
            file_data: UserFile = {
                "id": file_doc.id,
                "name": current_folder.name,
                "path": current_folder.parent.as_posix(),
                "size": 0,
                "extension": "folder",
                "isFolder": True,
                "storageLink": f"uploads/{user_id}{current_folder.as_posix()}",
                "uploadedOn": firestore.firestore.SERVER_TIMESTAMP,  # type: ignore
            }

            batch.create(file_doc, file_data)  # type: ignore
            count += 1
            # Firestore batches can only handle up to 500 operations
            if count == 500:
                batch.commit()
                batch = self.db.batch()
                count = 0

            # add parent to be created next
            folder_queue.append(current_folder.parent)

        if count > 0:
            batch.commit()

    def delete_path(self, user_id: str, path: str | Path) -> None:
        """
        Given a userID and path, deletes user file/folder matching the path.
        If the directory is a folder, deletes all files within the folder.

        Also deletes all corresponding files in Cloud Storage. (Empty folders
        may or not may not exist in Cloud Storage due to its flat structure.)
        """
        path = self.normalize_path(path)

        files_ref = self.db.collection("users", user_id, "files")
        bucket_ref = storage.Client().bucket("twocube-web.appspot.com")

        count = 0
        batch = self.db.batch()

        # find firestore doc corresponding to file/folder
        query = (
            files_ref.where(filter=FieldFilter("path", "==", path.parent.as_posix()))
            .where(filter=FieldFilter("name", "==", path.name))
            .limit(1)
        ).get()

        # if no doc found, we're done, unless we're deleting root
        if len(query) == 0 and path != Path("/"):
            return

        # if we're not deleting root, delete firestore record for that file
        is_folder = True
        if path != Path("/"):
            file_doc = query[0]
            is_folder: bool = file_doc.get("isFolder")
            storageLink: str = file_doc.get("storageLink")

            batch.delete(file_doc.reference)

            try:
                bucket_ref.blob(storageLink).delete()
            except NotFound:
                pass

        if is_folder:
            # delete firestore + cloud storage for all user files with given path
            # see here how i make this query: https://stackoverflow.com/a/46574143
            lex_next_string = path.as_posix()
            lex_next_string = lex_next_string[:-1] + chr(ord(lex_next_string[-1]) + 1)

            query = (
                files_ref.where(
                    filter=FieldFilter("path", ">=", path.as_posix())
                ).where(filter=FieldFilter("path", "<", lex_next_string))
            ).stream()

            for file_doc in query:
                storageLink = file_doc.get("storageLink")

                batch.delete(file_doc.reference)
                count += 1
                # Firestore batches can only handle up to 500 operations
                if count == 500:
                    batch.commit()
                    batch = self.db.batch()
                    count = 0

                try:
                    bucket_ref.blob(storageLink).delete()
                except NotFound:
                    pass

        if count > 0:
            batch.commit()
