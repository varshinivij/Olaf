import json

from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from ..datastructures.history import History
from ..services.agent_service import AgentService
from ..services.session_service import SessionService
from ..utils.validation import ValidationError, expect_values_in_request_body


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get"]))
def get_sessions(req: Request) -> Response:
    try:
        args = expect_values_in_request_body(req, ["userId"])
        user_id = args["userId"]

        session_service = SessionService()
        sessions = session_service.get_all_sessions(user_id)

        return Response(
            json.dumps(sessions),
            status=200,
            mimetype="application/json",
        )
    except ValidationError as e:
        return Response(
            json.dumps({"error": str(e)}), status=400, mimetype="application/json"
        )
    except Exception as e:
        return Response(
            json.dumps({"error": str(e)}), status=500, mimetype="application/json"
        )


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get"]))
def get_session_summary(req: Request) -> Response:
    try:
        args = expect_values_in_request_body(req, ["sessionId"])
        session_id = args["sessionId"]

        session_service = SessionService()
        session_data = session_service.get_session(session_id)

        if session_data is None:
            raise ValidationError("Session not found")

        agent_service = AgentService()
        history = History(session_data.get("history", []))

        summary_response = agent_service.chat_completion_summary(history)
        if summary_response is None:
            raise ValueError("Failed to generate summary")

        pdf_output = agent_service.create_pdf_from_summary(summary_response, session_id)

        return Response(
            pdf_output,
            content_type="application/pdf",
            status=200,
            headers={
                "Content-Disposition": f"attachment; filename={session_id}_summary.pdf"
            },
        )
    except ValidationError as e:
        return Response(
            json.dumps({"error": str(e)}), status=400, mimetype="application/json"
        )
    except Exception as e:
        return Response(
            json.dumps({"error": str(e)}), status=500, mimetype="application/json"
        )


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def rename_session(req: Request) -> Response:
    try:
        args = expect_values_in_request_body(req, ["sessionId", "newName"])
        session_id = args["sessionId"]
        new_name = args["newName"]

        session_service = SessionService()
        session_service.update_session(session_id, {"name": new_name})

        return Response(
            json.dumps({"message": "Session renamed successfully"}),
            status=200,
            mimetype="application/json",
        )
    except ValidationError as e:
        return Response(
            json.dumps({"error": str(e)}), status=400, mimetype="application/json"
        )
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["delete"]))
def delete_session(req: Request) -> Response:
    try:
        args = expect_values_in_request_body(req, ["sessionId"])
        session_id = args["sessionId"]

        session_service = SessionService()
        session_service.delete_session(session_id)

        return Response(
            json.dumps({"message": "Session deleted successfully"}),
            status=200,
            mimetype="application/json",
        )
    except ValidationError as e:
        return Response(
            json.dumps({"error": str(e)}), status=400, mimetype="application/json"
        )
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["delete"]))
def delete_all_sessions(req: Request) -> Response:
    try:
        args = expect_values_in_request_body(req, ["userId"])
        user_id = args["userId"]

        session_service = SessionService()
        session_service.delete_all_sessions(user_id)

        return Response(
            json.dumps({"message": "All sessions deleted successfully"}),
            status=200,
            mimetype="application/json",
        )
    except ValidationError as e:
        return Response(
            json.dumps({"error": str(e)}), status=400, mimetype="application/json"
        )
    except Exception as e:
        return Response(
            json.dumps({"error": str(e)}), status=500, mimetype="application/json"
        )
