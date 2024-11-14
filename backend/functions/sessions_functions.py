import firebase_admin
from firebase_admin import credentials, firestore
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
from functions_framework import http
import json
import flask

# Initialize Firebase app if not already initialized
if not firebase_admin._apps:
    cred = credentials.ApplicationDefault()
    firebase_admin.initialize_app(cred)

db = firestore.client()
    

# Delete Session
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["delete"]))
def delete_session(req: Request) -> Response:
    try:
        session_id = req.args.get('id')
        if not session_id:
            return Response(json.dumps({"error": "Session ID is required"}), status=400, mimetype='application/json')
        
        db.collection('sessions').document(session_id).delete()
        return Response(json.dumps({"message": "Session deleted successfully"}), status=200, mimetype='application/json')
    except Exception as e:
        print(f"Unable to delete session: {str(e)}")

# Get Sessions
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get"]))
def get_sessions(req: Request) -> Response:
    try:
        user_id = req.args.get('userId')
        if not user_id:
            return Response(json.dumps({"error": "User ID is required"}), status=400, mimetype='application/json')
        
        sessions_ref = db.collection('sessions').where('userId', '==', user_id)
        sessions = sessions_ref.stream()
        sessions_list = [{**doc.to_dict(), 'id': doc.id} for doc in sessions]
        return sessions_list
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')


# Delete All Sessions
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["delete"]))
def delete_all_sessions(req: Request) -> Response:
    try:
        user_id = req.args.get('userId')
        if not user_id:
            return Response(json.dumps({"error": "User ID is required"}), status=400, mimetype='application/json')
        
        sessions_ref = db.collection('sessions').where('userId', '==', user_id)
        sessions = sessions_ref.stream()
        
        batch = db.batch()
        count = 0
        for session in sessions:
            session_ref = db.collection('sessions').document(session.id)
            batch.delete(session_ref)
            count += 1
            # Firestore batches can only handle up to 500 operations
            if count == 500:
                batch.commit()
                batch = db.batch()
                count = 0
        if count > 0:
            batch.commit()
        
        return Response(json.dumps({"message": "All sessions deleted for user"}), status=200, mimetype='application/json')
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')
    

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def rename_session(req: Request) -> Response:
    try:
        session_id = req.args.get("session_id")
        new_name = req.json.get("newName")

        if not session_id:
            return Response(json.dumps({"error": "'session_id' is required"}), status=400)
        if not new_name:
            return Response(json.dumps({"error": "'newName' is required"}), status=400)

        # Update the session name
        db.collection("sessions").document(session_id).update({"name": new_name})
        return Response(json.dumps({"message": "Session renamed successfully"}), status=200, mimetype="application/json")

    except Exception as e:
        print(f"Error in rename_session: {str(e)}")
        return Response(json.dumps({"error": str(e)}), status=500)

def create_new_session(user_id, project_id, message, collection="sessions", ):
    """
    Create a new session in Firestore
    @return: session_id, session_data
    """
    new_session_ref = db.collection(collection).document()
    session_id = new_session_ref.id
    session_data = {
        "id": session_id,
        "userId": user_id,
        "projectId": project_id,
        "name": "<untitled session>",
        "context": "",
        "history": [
            {
                "type": "text",
                "role": "assistant",
                "content": "Hello, how can I help you today?"
            },
            {"type": "text", "role": "user", "content": message}
        ],
        "sandboxId": None
    }
    new_session_ref.set(session_data)
    return session_id,session_data

def update_session(session_id, session_data):
    """
    Update a session in Firestore
    @return: session_data
    """
    db.collection("sessions").document(session_id).update(session_data)
    session_ref = db.collection("sessions").document(session_id)
    session_data = session_ref.get().to_dict()
    return session_data

def add_message_to_session(session_id, message, role="user"):
    """
    Add a new message to the session history
    @return: session_data
    """
    db.collection("sessions").document(session_id).update({
        "history": firestore.ArrayUnion([{"type": "text", "role": role, "content": message}])
    })
    session_ref = db.collection("sessions").document(session_id)
    session_data = session_ref.get().to_dict()
    return session_data

def set_sand_box_id(session_id, sandbox_id):
    """
    Set the sandbox ID for a session
    @return: session_data
    """
    db.collection("sessions").document(session_id).update({
        "sandboxId": sandbox_id
    })
    session_ref = db.collection("sessions").document(session_id)
    session_data = session_ref.get().to_dict()
    return session_data
