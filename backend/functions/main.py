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
    try:
        sandbox_id = req.json.get("sandboxId")
        code = req.json.get("code")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        execution = sandbox.notebook.exec_cell(code)

        result = {"output": execution.logs.stdout}
        json_result = json.dumps(result)

        return Response(json_result, status=200)

    except Exception as e:
        json_error = json.dumps({"error": str(e)})
        return Response(json_error, status=500)
