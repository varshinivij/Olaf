import firebase_admin
from firebase_admin import credentials, firestore
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
from functions import Router
from functions.agent_utils import stream
from functions.agents.master_agent import MasterAgent
from functions.history import History
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
def master_agent_interaction(req: Request) -> Response:
    try:
        session_id = req.args.get("session_id")
        message = req.json.get("message")
        user_id = req.args.get('user_id')
        project_id = req.args.get('project_id')
        
        if not session_id:
            return flask.Response(json.dumps({"error": "'session_id' is required"}), status=400)
        if not message:
            return flask.Response(json.dumps({"error": "'message' is required"}), status=400)
        if not user_id:
            return flask.Response(json.dumps({"error": "'user_id' is required"}), status=400)
        if not project_id:
            return flask.Response(json.dumps({"error": "'project_id' is required"}), status=400)

        if session_id != "":
            session_data = db.collection("sessions").document(session_id).get().to_dict()
            if not session_data:
                return Response(json.dumps({"error": "Session not found"}), status=404, mimetype='application/json')
        else:
            new_session_ref = db.collection("sessions").document()
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
                    }
                ],
                "sandboxId": None
            }
            new_session_ref.set(session_data)
        
        session_data["history"].append({"role": "user", "content": message})

        router = Router()
        updated_session = router.route_session("Master", session_data)

        db.collection("sessions").document(session_id).set(updated_session, merge=True)

        assistant_message = updated_session["history"][-1]["content"]
        return Response(json.dumps({"reply": assistant_message}), status=200, mimetype='application/json')

    except Exception as e:
        print(f"Error in generate_plan: {str(e)}")
        return flask.Response(json.dumps({"error": str(e)}), status=500)