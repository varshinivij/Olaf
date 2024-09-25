import http
import time
from firebase_functions import https_fn, options
from e2b_code_interpreter import CodeInterpreter
from flask import send_file
from werkzeug.utils import secure_filename

from pathlib import Path
import io
import json

from google.cloud import storage
from flask import jsonify


E2B_API_KEY = "REMOVED"
E2B_TEMPLATE = "vh7kehbtf0t4xbx9ec9u"

@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["get"]))
def sandbox_status(sandbox_id: str) -> dict:
    """
    Returns the status of the E2B instance with the given sandbox ID.
    """
    try:
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        return {"alive": True,}
    except Exception as e:
        return {"alive": False, "error": str(e)}

@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def close_sandbox(req: https_fn.Request) -> https_fn.Response:
    """
    Closes the E2B instance with the given sandbox ID. This is useful for
    cleaning up resources and preventing memory leaks.

    Body data:
        sandboxId: the string ID of the sandbox to close
    """
    try:
        sandbox_id = req.json.get("sandboxId")
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        sandbox.close()

        return https_fn.Response("Sandbox closed", status=200)

    except Exception as e:
        json_error = json.dumps({"error": str(e)})
        return https_fn.Response(json_error, status=500)
    
    
@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))    
def request_sandbox(req: https_fn.Request) -> https_fn.Response:
    try:
        sandbox = CodeInterpreter(
            api_key=E2B_API_KEY,
            template=E2B_TEMPLATE,
            timeout=300,
            cwd="/home/user",
            # This is also the default folder for E2B uploads.
            # Setting the cwd might be a good idea since the root dir
            # contains the Dockerfile and whatnot. Maybe restrict user permits
            # to just this folder.
        )
        sandbox_id = sandbox.id
        print(sandbox_id)
        sandbox.keep_alive(5)  # keep box alive for 5minutes
        result = {"sandboxId": sandbox_id}
        json_result = json.dumps(result)

        return https_fn.Response(json_result, status=200)

    except Exception as e:
        return https_fn.Response({"error": str(e)}, status=500)


@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def upload_to_sandbox(req: https_fn.Request) -> https_fn.Response:
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
        sandbox_id = req.form["sandboxId"]
        files = req.files.getlist("file")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        for file in files:
            # Even though file is a FileStorage object (which is a custom
            # interface from Flask) and not a standard Python file object,
            # it still works. Ignore linters. This took me forever to figure out.
            if file.filename is not None:
                file.name = secure_filename(file.filename)
            # file.name will be "file" if no filename is provided
            # although i'm not sure if that's possible.

            sandbox.upload_file(file)

        # this is the directory where uploads are sent
        content = sandbox.filesystem.list("/home/user")
        cwdFiles = [{"name": item.name, "isDir": item.is_dir} for item in content]

        result = {"uploadDir": cwdFiles}

        json_result = json.dumps(result)
        return https_fn.Response(json_result, status=200)

    except Exception as e:
        json_error = json.dumps({"error": str(e)})
        return https_fn.Response(json_error, status=500)

@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def firebase_storage_to_sandbox(req: https_fn.Request) -> https_fn.Response:
    """
    Moves files from Firebase Storage to the E2B instance.
    """
    try:
        sandbox_id = req.json.get("sandboxId")
        file_paths = req.json.get("filePaths")

        if not sandbox_id or not file_paths:
            return jsonify({"error": "Missing sandboxId or filePaths"}), 400

        storage_client = storage.Client()
        bucket = storage_client.bucket("twocube-web.appspot.com")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)

        for file_path in file_paths:
            blob = bucket.blob(file_path)
            file_size = blob.size
            print(f"Downloading file: {file_path}, size: {file_size} bytes")
            file_name = file_path.split('/')[-1]
            blob.download_to_filename(file_name)

            #go from file name to file object
            file_name = Path(file_name)

            with open(file_name, "rb") as f:
                remote_path = sandbox.upload_file(f) 
                print(remote_path) 

        return jsonify({"message": "Files successfully moved to sandbox"}), 200
    
    except requests.exceptions.Timeout as e:
        # Handle timeout errors separately
        return jsonify({"error": "Request timed out", "details": str(e)}), 500

    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return https_fn.Response(json.dumps({"error": "An error occurred", "details": str(e)}), status=500)


@https_fn.on_request(cors=options.CorsOptions(cors_origins="*", cors_methods=["post"]))
def download_from_sandbox(req: https_fn.Request) -> https_fn.Response:
    """
    Downloads a file from an E2B instance relative to the "/home/user/" folder
    since that's where uploads go. It might be good to restrict permissions
    to just this folder since we don't want them messing with the Dockerfile
    and such.

    Because of this, I've set the CWD in the request_sandbox function to be
    "/home/user".

    Body data:
        sandboxId: the string ID of the sandbox to execute in
        path: the path inside the E2B instance to download. path is relative
              to "/home/user" (so you don't need to specify dir for uploads.)
    """
    try:
        sandbox_id = req.json.get("sandboxId")
        path: Path = Path("/home/user") / req.json.get("path")
        file_name = path.name

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        file_bytes = sandbox.download_file(path.as_posix())

        return send_file(
            io.BytesIO(file_bytes), as_attachment=True, download_name=file_name
        )

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

        print("Sandbox ID:", sandbox_id)

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
