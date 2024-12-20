from firebase_functions.https_fn import Request

def validate_master_agent_request(req: Request):
    session_id = req.args.get("session_id")
    message = req.args.get("message")
    user_id = req.args.get("user_id")
    project_id = req.args.get("project_id")

    if not message:
        raise ValueError("'message' is required")
    if not user_id:
        raise ValueError("'user_id' is required")
    if not project_id:
        raise ValueError("'project_id' is required")

    return session_id, user_id, project_id, message

def validate_name_maker_request(req: Request):
    # Expect JSON body with "history"
    data = req.json
    if not data or "history" not in data:
        raise ValueError("'history' field is required in request body")
    return data["history"]