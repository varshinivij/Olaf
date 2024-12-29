from typing import Any, Iterable

from firebase_admin import auth
from firebase_functions.https_fn import CallableRequest, Request
from werkzeug.datastructures import MultiDict


class ValidationError(Exception):
    """
    Validation error raised when a request fails validation.
    """

    def __init__(self, message: str):
        super().__init__(message)


def validate_request_id_token(req: Request) -> dict[str, Any]:
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
        raise ValidationError(
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
        raise ValidationError("'message' is required")
    if not user_id:
        raise ValidationError("'user_id' is required")
    if not project_id:
        raise ValidationError("'project_id' is required")

    return session_id, user_id, project_id, message


def validate_name_maker_request(req: Request):
    # Expect JSON body with "history"
    data = req.json
    if not data or "history" not in data:
        raise ValidationError("'history' field is required in request body")
    return data["history"]


def expect_values_in_request_args(
    req: Request, values: Iterable[str]
) -> MultiDict[str, str]:
    """
    Validates that all specified values are present in query args. Returns args
    after validating.
    """
    for value in values:
        if not req.args.get(value):
            raise ValidationError(f"'{value}' is required in request args")

    return req.args


def expect_values_in_request_body(
    req: Request, values: Iterable[str]
) -> dict[Any, Any]:
    """
    Validates that all specified values are present in request body as a JSON
    object. Returns body after validating.
    """
    data = req.json

    if not isinstance(data, dict):
        raise ValidationError("JSON object body expected")

    for value in values:
        if value not in data:
            raise ValidationError(f"'{value}' is required in request body")

    return data


def expect_values_in_request_data(
    req: CallableRequest, values: Iterable[str]
) -> dict[Any, Any]:
    """
    Validates that all specified values are present in request data.
    Returns body after validating. Only applicable to HTTP Callable functions.
    """
    data = dict(req.data)

    for value in values:
        if value not in data:
            raise ValidationError(f"'{value}' is required in request data")

    return data
