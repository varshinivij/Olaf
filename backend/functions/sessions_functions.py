import firebase_admin
from firebase_admin import credentials, firestore
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
from functions_framework import http
import json

# Initialize Firebase app if not already initialized
if not firebase_admin._apps:
    cred = credentials.ApplicationDefault()
    firebase_admin.initialize_app(cred)

db = firestore.client()

# Update Session
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["put"]))
def update_session(req: Request) -> Response:
    try:
        session_id = req.args.get('id')
        session_data = req.json
        if not session_id:
            return Response(json.dumps({"error": "Session ID is required"}), status=400, mimetype='application/json')
        
        db.collection('sessions2').document(session_id).set(session_data, merge=True)
        return Response(json.dumps({"message": "Session updated successfully"}), status=200, mimetype='application/json')
    except Exception as e:
        print(f"Unable to add session: {str(e)}")
    

# Delete Session
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["delete"]))
def delete_session(req: Request) -> Response:
    try:
        session_id = req.args.get('id')
        if not session_id:
            return Response(json.dumps({"error": "Session ID is required"}), status=400, mimetype='application/json')
        
        db.collection('sessions2').document(session_id).delete()
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
        
        sessions_ref = db.collection('sessions2').where('userId', '==', user_id)
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
        
        sessions_ref = db.collection('sessions2').where('userId', '==', user_id)
        sessions = sessions_ref.stream()
        
        batch = db.batch()
        count = 0
        for session in sessions:
            session_ref = db.collection('sessions2').document(session.id)
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
def create_session(req: Request) -> Response:
    try:
        session = req.json
        doc_ref = db.collection('sessions2').add(session)
        return Response(json.dumps({"id": doc_ref[1].id}), status=200, mimetype='application/json')
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')
