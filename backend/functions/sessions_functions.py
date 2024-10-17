from firebase_admin import firestore
from firebase_functions.https_fn import Request, Response
import json

db = firestore.client()

def add_session(session) -> None:
    try:
        db.collection('sessions2').add(session)
    except Exception as e:
        print(f"Unable to add session: {str(e)}")

def update_session(session_id, session_data) -> None:
    try:
        db.collection('sessions2').document(session_id).set(session_data, merge=True)
    except Exception as e:
        print(f"Unable to update session: {str(e)}")

def delete_session(session_id) -> None:
    try:
        db.collection('sessions2').document(session_id).delete()
    except Exception as e:
        print(f"Unable to delete session: {str(e)}")

def get_sessions(user_id) -> list:
    try:
        sessions_ref = db.collection('sessions2').where('userId', '==', user_id)
        sessions = sessions_ref.stream()
        sessions_list = [{**doc.to_dict(), 'id': doc.id} for doc in sessions]
        return sessions_list
    except Exception as e:
        print(f"Unable to get sessions: {str(e)}")
        sessions_list = []
        
