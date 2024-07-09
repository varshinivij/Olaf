from firebase_functions import https_fn, options
from firebase_admin import initialize_app
from agent_logic import MasterAgent
from e2b import Sandbox
import json

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
        sandbox = Sandbox(api_key=e2b_api_key, template="base", timeout=300)
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

        sandbox = Sandbox.reconnect(sandbox_id, api_key=e2b_api_key) 
        echo = sandbox.process.start('echo "Code Ran (not really, this is a place holder)!"')  
        echo.wait()
        result = {"output": echo.stdout}
        json_result = json.dumps(result)
        return https_fn.Response(json_result, status=200)
    except Exception as e:
        json_error = json.dumps({'error': str(e)})
        return https_fn.Response(json_error, status=500)

