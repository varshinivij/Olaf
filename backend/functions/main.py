from firebase_functions import https_fn, storage_fn, firestore_fn, options
from firebase_admin import initialize_app, firestore
from agent_logic import MasterAgent
from e2b import Sandbox
from e2b_code_interpreter import CodeInterpreter
import json
import re

initialize_app()

# --- EXAMPLE ---
@https_fn.on_request(
    cors=options.CorsOptions(cors_origins="*", cors_methods=["get", "post"])
)
def on_request_example(req: https_fn.Request) -> https_fn.Response:
    return https_fn.Response("Hello world this is me!!")

# --- LLM Agent Function ---
@https_fn.on_request(
    cors=options.CorsOptions(cors_origins="*", cors_methods=["post"])
)
def ask_agent(req: https_fn.Request) -> https_fn.Response:
    #a request will be conversation history
    history = req.json.get("history")
    agent = MasterAgent()
    response = agent.chat_completion(history)
    parsed_response = agent.parse_output(response)
    return https_fn.Response(parsed_response)

# --- E2B Function ---
e2b_api_key = 'REMOVED'

@https_fn.on_request(
    cors=options.CorsOptions(cors_origins="*", cors_methods=["post"])
)
def request_sandbox(req: https_fn.Request) -> https_fn.Response:
    try:
        sandbox = CodeInterpreter(api_key=e2b_api_key, template="vh7kehbtf0t4xbx9ec9u", timeout=300)
        sandbox_id = sandbox.id
        sandbox.keep_alive(5 * 60) #keep box alive for 5minutes
        result = {"sandboxId": sandbox_id}
        json_result = json.dumps(result)
        return https_fn.Response(json_result, status=200)
    except Exception as e:
        return https_fn.Response({'error': str(e)}, status=500)

@https_fn.on_request(
    cors=options.CorsOptions(cors_origins="*", cors_methods=["post"])
)
def execute_on_sandbox(req: https_fn.Request) -> https_fn.Response:
    try:
        sandbox_id = req.json.get("sandboxId")
        code = req.json.get("code")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=e2b_api_key)
        execution = sandbox.notebook.exec_cell(code)
        result = {"output": execution.logs.stdout}
        json_result = json.dumps(result)
        return https_fn.Response(json_result, status=200)
    except Exception as e:
        json_error = json.dumps({'error': str(e)})
        return https_fn.Response(json_error, status=500)


# --- Uploading Files to Cloud Storage Updates user/ID/fields ---
@storage_fn.on_object_finalized()
def on_file_uploaded(event: storage_fn.CloudEvent) -> None:
    file_path = event.data.name

    match = re.match(r'^uploads/([^/]+)/(.+)$', file_path)
    if not match:
        return

    user_uid, file_name = match.groups()

    db = firestore.client()
    user_ref = db.collection('users').document(user_uid)

    # Get the current user document
    user_doc = user_ref.get()

    if not user_doc.exists:
        print(f"User document for user {user_uid} does not exist.")
        return
    else:
        current_files = user_doc.to_dict().get('files', [])
        if not isinstance(current_files, list):
            print(f"'files' field for user {user_uid} is not an array. Overwriting with new array.")
            user_ref.update({'files': []})

        if file_path not in current_files:
            # Update the user's document
            user_ref.update({
                'files': firestore.ArrayUnion([file_path])
            })

            print(f"Added {file_path} to files array for user {user_uid}")
