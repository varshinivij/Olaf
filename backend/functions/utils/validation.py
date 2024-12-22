from firebase_admin import auth
from firebase_functions.https_fn import Request


def validate_request_id_token(req: Request) -> dict:
    """
    Validates the Firebase ID token passed in a request's Authorization header.
    Returns the decoded token as a dict containing user data. Example returned
    payload:
    {
        "name": str
        "picture": str (link)
        "iss": str (link),
        "aud": str,
        "auth_time": int,
        "user_id": str,
        "sub": str,
        "iat": int,
        "exp": int,
        "email": str,
        "email_verified": bool,
        "firebase": {
            "identities": {
                "google.com": list[str],
                "email": list[str],
            },
            "sign_in_provider": str,
        },
        "uid": str,
    }

    To send the ID token in Angular:
    const headers = new HttpHeaders({
      Authorization: `Bearer ${await this.auth.currentUser?.getIdToken()}`,
    });
    this.http.post(url, { param1: 'folder1' }, { headers });

    Not necessary for HTTP Callables, only pure HTTP functions.
    """
    if not req.headers.get("Authorization", "").startswith("Bearer "):
        raise PermissionError(
            "No Firebase ID token passed as Bearer token in Authorization header."
            " Provide the following HTTP header: Authorization: Bearer <Firebase ID Token>"
        )

    id_token = req.headers["Authorization"].split("Bearer ")[1]
    return auth.verify_id_token(id_token, check_revoked=True)


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
