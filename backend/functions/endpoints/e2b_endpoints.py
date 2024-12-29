import io
import json
import requests

from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
from flask import send_file

from functions.services.e2b_service import E2BService
from functions.utils.validation import ValidationError, expect_values_in_request_body


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get"]))
def sandbox_status(req: Request) -> Response:
    """
    Returns the status of the E2B instance with the given sandbox ID.
    """
    try:
        args = expect_values_in_request_body(req, ["sandboxId"])
        sandbox_id = args["sandboxId"]

        e2b_service = E2BService()
        alive = e2b_service.get_is_sandbox_alive(sandbox_id)

        return Response(
            json.dumps({"alive": alive}),
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


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def close_sandbox(req: Request) -> Response:
    """
    Closes the E2B instance with the given sandbox ID. This is useful for
    cleaning up resources and preventing memory leaks.
    """
    try:
        args = expect_values_in_request_body(req, ["sandboxId"])
        sandbox_id = args["sandboxId"]

        e2b_service = E2BService()
        e2b_service.close_sandbox(sandbox_id)

        return Response(
            json.dumps({"message": "Sandbox closed successfully"}),
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


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def request_sandbox(req: Request) -> Response:
    try:
        e2b_service = E2BService()
        sandbox_id = e2b_service.create_sandbox()

        return Response(
            json.dumps({"sandboxId": sandbox_id}),
            status=200,
            mimetype="application/json",
        )
    except Exception as e:
        return Response(
            json.dumps({"error": str(e)}), status=500, mimetype="application/json"
        )


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def upload_to_sandbox(req: Request) -> Response:
    """
    Accepts HTTP requests with ContentType "multipart/form-data." This can be
    done through <form> tags with enctype="multipart/form-data" or using
    FormData objects: https://developer.mozilla.org/en-US/docs/Web/API/FormData

    Currently, the form data should have a key "sandboxId": the key to the
    sandbox, along with any files that should be uploaded to the E2B instance
    with FormData key as "file".

    For future consideration, it may be better to use a RESTful interface:
    for example, [api]/sandbox/{sandboxId}/upload. But for now, the ID goes
    inside the form data.
    """
    try:
        if "sandboxId" not in req.form:
            raise ValidationError("Missing sandboxId in form data")

        sandbox_id = req.form["sandboxId"]
        files = req.files.getlist("file")

        e2b_service = E2BService()
        e2b_service.upload_formdata_files_to_sandbox(sandbox_id, files)

        return Response(
            json.dumps({"message": "Files uploaded successfully"}),
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


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def firebase_storage_to_sandbox(req: Request) -> Response:
    """
    Moves files from Firebase Storage to the E2B instance.
    """
    try:
        args = expect_values_in_request_body(req, ["sandboxId", "filePaths"])
        sandbox_id = args["sandboxId"]
        file_paths = args["filePaths"]

        e2b_service = E2BService()
        e2b_service.upload_firebase_files_to_sandbox(sandbox_id, file_paths)

        return Response(
            json.dumps({"message": "Files uploaded successfully"}),
            status=200,
            mimetype="application/json",
        )
    except requests.exceptions.Timeout as e:
        # Handle timeout errors separately
        return Response(
            {"error": "Request timed out", "details": str(e)},
            status=500,
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


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def download_from_sandbox(req: Request) -> Response:
    """
    Downloads a file from an E2B instance relative to the "/home/user/" folder
    since that's where uploads go. It might be good to restrict permissions
    to this folder since we don't want them messing with the Dockerfile, etc.

    Because of this, I've set the CWD in the request_sandbox function to be
    "/home/user".

    Body data:
        sandboxId: the string ID of the sandbox to execute in
        path: the path inside the E2B instance to download. path is relative
              to "/home/user" (so you don't need to specify dir for uploads.)
    """
    try:
        args = expect_values_in_request_body(req, ["sandboxId", "path"])
        sandbox_id = args["sandboxId"]
        path = args["path"]

        e2b_service = E2BService()
        file_name, file_bytes = e2b_service.download_file_from_sandbox(sandbox_id, path)

        return send_file(
            io.BytesIO(file_bytes), as_attachment=True, download_name=file_name
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
def execute_on_sandbox(req: Request) -> Response:
    """
    This endpoint essentially returns the properties of the E2B Execution object
    after executing some code.

    Body data:
        sandboxId: the string ID of the sandbox to execute in
        code: the Python code to execute as a string

    Returns:
        logs: {
            stdout: List[str] of stdout messages
            stderr: List[str] of stderr messages
        }
        error: None or str of E2B Error message
        results: List[ {mimetype: str of data or None} ]
            example
            List [{
                    "image/png": (base64 encoded PNG data as string)
                    "text/plain": "<Figure size 1000x600 with 1 Axes>"
                },
                ...]

    According to the E2B docs https://e2b.dev/docs/hello-world/py
    the image data is encoded as b64, which can then be used as a binary string.
    Text is returned as plain text. Not sure about other files like PDFs, docs
    don't specify.

    Use the image string directly in an <img> tag's src attribute and specify
    it as b64 to display it. The conversion is built into HTML already.

    I don't see easier ways to return multiple files at once besides zipping,
    but then I don't think you can send additional JSON data along with it.
    Raw binary strings from E2B are the most convenient.

    A single result could have multiple mimetypes. I will return as much data
    as possible without filtering since it's the frontend's job to filter data
    and the backend's to give as much as possible.
    """
    try:
        args = expect_values_in_request_body(req, ["sandboxId", "code"])
        sandbox_id = args["sandboxId"]
        code = args["code"]

        e2b_service = E2BService()
        result = e2b_service.execute_code_in_sandbox(sandbox_id, code)

        return Response(
            json.dumps(result),
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


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def run_terminal_command(req: Request) -> Response:
    """
    Runs a terminal command on the E2B sandbox terminal using exec_cell.

    Body data:
        sandboxId: the string ID of the sandbox to connect to
        command: the terminal command to execute

    Returns:
        stdout: The standard output of the command
        stderr: The standard error of the command
        error: Any traceback or error messages (if applicable)
    """
    try:
        args = expect_values_in_request_body(req, ["sandboxId", "command"])
        sandbox_id = args["sandboxId"]
        command = args["command"]

        e2b_service = E2BService()
        result = e2b_service.execute_command_in_sandbox(sandbox_id, command)

        return Response(
            json.dumps(result),
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
