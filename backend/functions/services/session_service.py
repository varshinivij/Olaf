from typing import Tuple, Dict, Any
from firebase_admin import firestore
import sessions_functions

db = firestore.client()

def get_or_create_session(user_id: str, project_id: str, message: str, session_id: str = None) -> Tuple[str, Dict[str, Any]]: # type: ignore
    """
    If session_id is provided, append the user message. Otherwise, create a new session.
    Returns the session_id and session_data dictionary.
    """
    if session_id:
        session_data = sessions_functions.add_message_to_session(session_id, message, role="user")
    else:
        session_id, session_data = sessions_functions.create_new_session(user_id, project_id, message)
    if not session_data:
        raise ValueError("Session not found or could not be created.")
    return session_id, session_data

def append_message_to_session(session_id: str, message: str, role: str = "assistant") -> None:
    sessions_functions.add_message_to_session(session_id, message, role=role)