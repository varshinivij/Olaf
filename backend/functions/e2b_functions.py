from firebase_functions import https_fn, options
from e2b_code_interpreter import CodeInterpreter

import json


E2B_API_KEY = "REMOVED"
E2B_TEMPLATE = "vh7kehbtf0t4xbx9ec9u"


@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def request_sandbox(req: https_fn.Request) -> https_fn.Response:
    try:
        sandbox = CodeInterpreter(
            api_key=E2B_API_KEY, template=E2B_TEMPLATE, timeout=300
        )
        sandbox_id = sandbox.id
        sandbox.keep_alive(5 * 60)  # keep box alive for 5minutes
        result = {"sandboxId": sandbox_id}
        json_result = json.dumps(result)

        return https_fn.Response(json_result, status=200)

    except Exception as e:
        return https_fn.Response({"error": str(e)}, status=500)


@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def upload_to_sandbox(req: https_fn.Request) -> https_fn.Response:
    """
    Body data:
        sandboxId: the string ID of the sandbox to execute in

    Uploaded files should come from a <form enctype="multipart/form-data">
    with input tags for file upload. That seems to be how Flask works, not my
    design. (Firebase Cloud Functions are built on top of Flask.) View the
    documentation under req.files (Intellisense helps) to learn more.
    """
    try:
        sandbox_id = req.json.get("sandboxId")
        files = req.files.getlist("files")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        for file in files:
            sandbox.upload_file(file.stream)

        content = sandbox.filesystem.list(".")
        cwdFiles = [{"name": item.name, "isDir": item.is_dir} for item in content]

        result = {"cwd": cwdFiles}

        json_result = json.dumps(result)
        return https_fn.Response(json_result, status=200)

    except Exception as e:
        json_error = json.dumps({"error": str(e)})
        return https_fn.Response(json_error, status=500)


@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def execute_on_sandbox(req: https_fn.Request) -> https_fn.Response:
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
            example:
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

    # to put auth into this endpoint, use the validate_firebase_id_token()
    # function in file_storage_functions and wrap it around the req: Request
    # object. this will require that the request puts a
    # "Authorization: Bearer {idToken}" with an idToken from Firebase or it
    # will error. Not currently putting it since it might be better to put
    # the function somewhere else and don't want to login in everytime while
    # testing.

    try:
        sandbox_id = req.json.get("sandboxId")
        code = req.json.get("code")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        execution = sandbox.notebook.exec_cell(code)

        result = {
            "logs": {"stdout": execution.logs.stdout, "stderr": execution.logs.stderr},
            "error": execution.error.traceback if execution.error else None,
            "results": [result.raw for result in execution.results],
        }

        json_result = json.dumps(result)

        return https_fn.Response(json_result, status=200)

    except Exception as e:
        json_error = json.dumps({"error": str(e)})
        return https_fn.Response(json_error, status=500)
