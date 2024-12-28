from typing import List, Optional

from firebase_admin import firestore
from google.cloud.firestore import ArrayUnion

from functions.models.chat_message import ChatMessage, ChatMessageRole
from functions.models.session import Session


class SessionService:
    """
    Service class containing utility methods for session endpoints.
    """

    def __init__(self):
        self.db = firestore.client()

    def add_message_to_session(
        self, session_id: str, message: str, role: ChatMessageRole
    ) -> Session | None:
        """
        Adds a message to a given session in Firestore.
        """
        return self.update_session(
            session_id,
            {
                "history": ArrayUnion(
                    [{"type": "text", "role": role, "content": message}]
                )
            },
        )

    def create_session(
        self, user_id: str, project_id: str, message: Optional[ChatMessage] = None
    ) -> Session:
        """
        Creates a new session in Firestore. Optionally add a message upon creation.
        """
        new_session_ref = self.db.collection("sessions").document()

        session_data: Session = {
            "id": new_session_ref.id,
            "userId": user_id,
            "projectId": project_id,
            "name": "<untitled session>",
            "context": "",
            "history": [
                {
                    "type": "text",
                    "role": "assistant",
                    "content": "Hello, how can I help you today?",
                },
            ],
            "sandboxId": None,
        }

        if message is not None:
            session_data["history"].append(message)

        new_session_ref.set(session_data)  # type: ignore
        return session_data

    def get_session(self, session_id: str) -> Session | None:
        """
        Gets a given session from Firestore.
        """
        session_ref = self.db.collection("sessions").document(session_id)
        session_data = session_ref.get().to_dict()
        return session_data  # type: ignore

    def get_all_sessions(self, user_id: str) -> List[Session]:
        """
        Gets all sessions for a given user from Firestore.
        """
        sessions_ref = self.db.collection("sessions").where("userId", "==", user_id)
        sessions_list: List[Session] = []
        sessions = sessions_ref.stream()

        for doc in sessions:
            session_data = doc.to_dict()
            session_data["id"] = doc.id
            sessions_list.append(session_data)

        return sessions_list

    def update_session(self, session_id: str, session_data: dict) -> Session | None:
        """
        Updates a session in Firestore.
        """
        self.db.collection("sessions").document(session_id).update(session_data)
        session_ref = self.db.collection("sessions").document(session_id)
        session_data = session_ref.get().to_dict()  # type: ignore
        return session_data  # type: ignore

    def delete_session(self, session_id: str) -> None:
        """
        Deletes a session in Firestore.
        """
        self.db.collection("sessions").document(session_id).delete()

    def delete_all_sessions(self, user_id: str) -> None:
        """
        Deletes all sessions for a given user in Firestore.
        """
        sessions_ref = self.db.collection("sessions").where("userId", "==", user_id)
        sessions = sessions_ref.stream()

        count = 0
        batch = self.db.batch()
        for doc in sessions:
            session_ref = self.db.collection("sessions").document(doc.id)
            batch.delete(session_ref)
            count += 1
            # Firestore batches can only handle up to 500 operations
            if count == 500:
                batch.commit()
                batch = self.db.batch()
                count = 0

        if count > 0:
            batch.commit()
