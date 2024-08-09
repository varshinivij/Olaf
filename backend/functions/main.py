from firebase_admin import initialize_app

from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from e2b import Sandbox
from e2b_code_interpreter import CodeInterpreter

from agent_logic import MasterAgent

import json

# initialize cloud functions defined in this file  (organizing like this
# may be better). maybe even put helper modules in a separate /utils folder.
from file_storage_functions import (
    handle_user_file_upload,
    request_user_create_folder,
    request_user_delete_path,
)


initialize_app()


E2B_API_KEY = "REMOVED"


# --- EXAMPLE ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get", "post"]))
def on_request_example(req: Request) -> Response:
    return Response("Hello world this is me!!")


# --- LLM Agent Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def ask_agent(req: Request) -> Response:
    # a request will be conversation history
    history = req.json.get("history")
    agent = MasterAgent()
    response = agent.chat_completion(history)
    parsed_response = agent.parse_output(response)

    return Response(parsed_response)


# --- E2B Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def request_sandbox(req: Request) -> Response:
    try:
        sandbox = CodeInterpreter(
            api_key=E2B_API_KEY, template="vh7kehbtf0t4xbx9ec9u", timeout=300
        )
        sandbox_id = sandbox.id
        sandbox.keep_alive(5 * 60)  # keep box alive for 5minutes
        result = {"sandboxId": sandbox_id}
        json_result = json.dumps(result)

        return Response(json_result, status=200)

    except Exception as e:
        return Response({"error": str(e)}, status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def execute_on_sandbox(req: Request) -> Response:
    """
    This endpoint essentially returns the properties of the E2B Execution object
    after executing some code.

    Returns:
        logs: {
            stdout: List[str] of stdout messages
            stderr: List[str] of stderr messages
        }
        error: None or str of E2B Error message
        results: List[ {mimetype: str of data or None} ]
            example:
            List [
                {
                    "image/png": (base64 encoded PNG data as string)
                }
            ]

    According to the E2B docs https://e2b.dev/docs/hello-world/py
    the image data is encoded as b64, which can then be used as a binary string.
    Not currently sure whether all data is encoded as b64 since the docs are
    unclear. Will experiment to find out.

    Decode the strings from B64 on the frontend and make the image from there
    as a binary string.

    I don't see easier ways to return multiple files at once besides zipping,
    but then I don't think you can send additional JSON data along with it.
    Raw binary strings strings from E2B are the most convenient.

    Also unclear whether a single result could have multiple mimetypes, but I
    will return as much data as possible without filtering since it's the
    frontend's job to filter data and the backend's to give as much as possible.
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

        return Response(json_result, status=200)

    except Exception as e:
        json_error = json.dumps({"error": str(e)})
        return Response(json_error, status=500)
