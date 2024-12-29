from flask import stream_with_context
import json

from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from ..datastructures.history import History
from ..datastructures.router import Router
from ..routes.masteragent_route import masteragent_route
from ..services.agent_service import chat_completion
from ..services.session_service import SessionService
from ..utils.validation import (
    validate_master_agent_request,
    validate_name_maker_request,
)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def master_agent_interaction(req: Request) -> Response:
    try:
        session_service = SessionService()

        # Extract necessary parameters and validate them
        session_id, user_id, project_id, message = validate_master_agent_request(req)
        if not session_id:
            session_id = None

        # Retrieve or create session with the new user message
        if session_id is None:
            session_data = session_service.create_session(
                user_id,
                project_id,
                {"type": "text", "role": "user", "content": message},
            )
            session_id = session_data["id"]
        else:
            session_data = session_service.add_message_to_session(
                session_id, {"type": "text", "role": "user", "content": message}
            )

        if session_data is None:
            raise ValueError("Session not found")

        # Filter out non-text/code/error messages
        # NOTE this logic can be put in History class
        session_data["history"] = [
            m
            for m in session_data["history"]
            if m.get("type") in ["text", "code", "error", "result", "executedCode"]
        ]

        router = Router()
        router.add_route("master", masteragent_route)
        response_generator = router.route("master", session_data["history"])  # type: ignore

        # NOTE can probably be put into agent service.
        # the add_message_to_session() call can probably be done by putting a
        # callback parameter in AgentService.generate_stream(). this would
        # let you keep the SessionService.add_message_to_session() call here
        # without having to put a SessionService object to the agent_service.

        def generate_stream():
            full_response = ""
            for chunk in response_generator:
                try:
                    data = json.dumps(
                        {"type": chunk["type"], "content": chunk["content"]}
                    )
                    yield f"data: {data}\n\n".encode("utf-8")
                    full_response += chunk["content"]
                except Exception as e:
                    # Log and skip malformed chunk
                    print(f"Error streaming chunk: {e}")
                    continue
            # Store the assistant's response
            session_service.add_message_to_session(
                session_id,
                {"type": "text", "role": "assistant", "content": full_response},
            )

        return Response(
            stream_with_context(generate_stream()),
            headers={"session_id": session_id},
            mimetype="text/event-stream",
        )

    except Exception as e:
        print(f"Error in master_agent_interaction_handler: {str(e)}")
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def name_maker(req: Request) -> Response:
    try:
        history = validate_name_maker_request(req)

        system = """
        You are given a conversation history.
        Please construct a super short title for the conversation.
        """
        hist = History(history)
        hist.log("system", system, "text")

        response = chat_completion(hist)
        response_str = ""
        for chunk in response:
            response_str += (
                chunk.get("choices", [{}])[0].get("delta", {}).get("content", "")
            )

        response_data = {"message": response_str}
        return Response(json.dumps(response_data), status=200)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)
