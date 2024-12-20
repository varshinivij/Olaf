import json
from flask import Response, stream_with_context
from firebase_functions.https_fn import Request
from services.session_service import get_or_create_session, append_message_to_session
from services.router_service import get_master_router
from utils.validation import validate_master_agent_request

def master_agent_interaction_handler(req: Request) -> Response:
    try:
        # Extract necessary parameters and validate them
        session_id, user_id, project_id, message = validate_master_agent_request(req)
        if not session_id:
            session_id = None
        # Retrieve or create session with the new user message
        session_id, session_data = get_or_create_session(user_id, project_id, message, session_id) # type: ignore

        # Filter out non-text/code/error messages
        session_data["history"] = [
            m for m in session_data["history"] 
            if m.get("type") in ["text", "code", "error", "result", "executedCode"]
        ]

        router = get_master_router()
        response_generator = router.route("master", session_data["history"]) # type: ignore

        def generate_stream():
            full_response = ""
            for chunk in response_generator:
                try:
                    data = json.dumps({"type": chunk["type"], "content": chunk["content"]})
                    yield f"data: {data}\n\n".encode('utf-8')
                    full_response += chunk["content"]
                except Exception as e:
                    # Log and skip malformed chunk
                    print(f"Error streaming chunk: {e}")
                    continue
            # Store the assistant's response
            append_message_to_session(session_id, full_response, role="assistant")

        return Response(stream_with_context(generate_stream()), headers={"session_id": session_id}, mimetype="text/event-stream")

    except Exception as e:
        print(f"Error in master_agent_interaction_handler: {str(e)}")
        return Response(json.dumps({"error": str(e)}), status=500)