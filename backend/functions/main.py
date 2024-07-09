from firebase_functions import https_fn, options
from firebase_admin import initialize_app
from agent_logic import MasterAgent

initialize_app()

@https_fn.on_request(
    cors=options.CorsOptions(cors_origins="*", cors_methods=["get", "post"])
)
def on_request_example(req: https_fn.Request) -> https_fn.Response:
    return https_fn.Response("Hello world this is me!!")


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

