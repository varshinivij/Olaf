from firebase_admin import initialize_app

from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from agent_logic import MasterAgent


initialize_app()

# import other modules' Cloud Functions
from e2b_functions import (
    request_sandbox,
    execute_on_sandbox,
    upload_to_sandbox,
    download_from_sandbox
)

from file_storage_functions import (
    handle_user_file_upload,
    request_user_create_folder,
    request_user_delete_path,
)


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
