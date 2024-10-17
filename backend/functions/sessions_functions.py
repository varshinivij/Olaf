import firebase_admin
from firebase_admin import credentials, firestore
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
from firebase_admin import initialize_app
from functions_framework import http
import json

db = firestore.client()

def add_session(req: Request) -> Response:
    try:
        session = req.json
        doc_ref = db.collection('sessions2').add(session)
        return Response(json.dumps({"id": doc_ref[1].id}), status=200, mimetype='application/json')
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')

def update_session(req: Request) -> Response:
    try:
        session_id = req.args.get('id')
        session_data = req.json
        if not session_id:
            return Response(json.dumps({"error": "Session ID is required"}), status=400, mimetype='application/json')
        
        db.collection('sessions2').document(session_id).set(session_data, merge=True)
        return Response(json.dumps({"message": "Session updated successfully"}), status=200, mimetype='application/json')
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')

def delete_session(req: Request) -> Response:
    try:
        session_id = req.args.get('id')
        if not session_id:
            return Response(json.dumps({"error": "Session ID is required"}), status=400, mimetype='application/json')
        
        db.collection('sessions2').document(session_id).delete()
        return Response(json.dumps({"message": "Session deleted successfully"}), status=200, mimetype='application/json')
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')

def get_sessions(req: Request) -> Response:
    try:
        user_id = req.args.get('userId')
        if not user_id:
            return Response(json.dumps({"error": "User ID is required"}), status=400, mimetype='application/json')
        
        sessions_ref = db.collection('sessions2').where('userId', '==', user_id)
        sessions = sessions_ref.stream()
        sessions_list = [{**doc.to_dict(), 'id': doc.id} for doc in sessions]
        
        return Response(json.dumps(sessions_list), status=200, mimetype='application/json')
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500, mimetype='application/json')
