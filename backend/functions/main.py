import json
import time
from typing import Callable

from executor import Executor
from agents.codemaster_agent import CodeMasterAgent
import flask

from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
import concurrent.futures
from history import History
from functions_framework import http
from agent_utils import chat_completion, stream
from firebase_admin import firestore
initialize_app()
import sessions_functions
from router import Router
from pipe import Pipe

# import other modules' Cloud Functions
from e2b_functions import (
    request_sandbox,
    execute_on_sandbox,
    upload_to_sandbox,
    download_from_sandbox,
    sandbox_status,
    close_sandbox,
    firebase_storage_to_sandbox
)

from file_storage_functions import (
    handle_user_file_upload,
    request_user_create_folder,
    request_user_delete_path
)

from sessions_functions import (
    add_message_to_session,
    delete_session,
    get_sessions,
    delete_all_sessions,
    rename_session,
    get_session_summary
)

# This needs to be cleaned up
E2B_API_KEY = "REMOVED"
db = firestore.client()

def master_route_function(session):
    history = session
    if not history:
            return None
    history = History(history)
    master_agent = CodeMasterAgent("python",history) # type: ignore #using CodeMasterAgent
    return master_agent.generate_response()
    

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def master_agent_interaction(req: Request) -> Response:
    try:
        session_id = req.args.get("session_id")
        message = req.args.get("message")
        user_id = req.args.get("user_id")
        project_id = req.args.get("project_id")
        print(f"session_id: {session_id}")
        print(f"message: {message}")
        print(f"user_id: {user_id}")
        print(f"project_id: {project_id}")
        if not message:
            return Response(json.dumps({"error": "'message' is required"}), status=400)
        if not user_id:
            return Response(json.dumps({"error": "'user_id' is required"}), status=400)
        if not project_id:
            return Response(json.dumps({"error": "'project_id' is required"}), status=400)
        if session_id: 
            # Append the new user message to the existing session
            print("adding message to session")
            session_data = sessions_functions.add_message_to_session(session_id, message, role="user")
            print("added message to session")
        else:
            session_id, session_data = sessions_functions.create_new_session(user_id, project_id, message)
        
        if not session_data:
            return Response(json.dumps({"error": "Session not found"}), status=404)
        # Filter session so that history is only made of text, code, and errors (no images)
        session_data["history"] = [message for message in session_data["history"] if message["type"] in ["text", "code", "error", "result", "executedCode"]]
        # Route the session history for response generation
        router = Router()
        print("created router")
        router.add_route("master", master_route_function)
        response_generator = router.route("master", session_data["history"]) # type: ignore
        print(f"response_generator: {response_generator}")
        def generate_stream():
            full_response = ""
            for chunk in response_generator:
                try:
                    # Prepare a JSON object with 'type' and 'content'
                    data = json.dumps({
                        "type": chunk["type"],     # 'text' or 'code'
                        "content": chunk["content"]
                    })
                    yield f"data: {data}\n\n".encode('utf-8')
                    full_response += chunk["content"]
                except Exception as e:
                    # Handle malformed chunk gracefully
                    print(f"Malformed chunk: {chunk} - Error: {str(e)}")
                    continue
            # Update the session history with the assistant's full response
            sessions_functions.add_message_to_session(session_id, full_response, role="assistant")
        return Response(flask.stream_with_context(generate_stream()), headers={"session_id": session_id}, mimetype="text/event-stream")   
    except Exception as e:
        print(f"Error in master_agent_interaction: {str(e)}")
        return flask.Response(json.dumps({"error": str(e)}), status=500)

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def name_maker(req: Request) -> Response:
    # a request will be conversation history
    system = """
    You are given a conversation history.
    Please construct a super short title for the conversation.
    """
    history = req.json.get("history")
    history = History(history)
    history.log("system", system)
    response = chat_completion(history)
    # Initialize an empty string to accumulate the response text
    response_str = ""
    for chunk in response:
        response_str += chunk.get("choices", [{}])[0].get("delta", {}).get("content", "")
    response_data = {
        "message": response_str,
    }
    # Return the response as a JSON response
    return Response(json.dumps(response_data), status=200)

