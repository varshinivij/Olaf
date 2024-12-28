from firebase_admin import firestore


class FileStorageService:
    """
    Service class containing utility methods for file storage endpoints.

    Currently all of the logic remains in file_storage_endpoints.py.
    """

    def __init__(self):
        self.db = firestore.client()
