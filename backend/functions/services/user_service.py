from firebase_admin import firestore


class UserService:
    """
    Service class containing utility methods for user endpoints.
    """

    def __init__(self):
        self.db = firestore.client()

    def get_user_exists(self, user_id: str) -> bool:
        """
        Returns whether a user of the specified ID exists.
        """
        user_doc = self.db.collection("users").document(user_id).get()
        return user_doc.exists
