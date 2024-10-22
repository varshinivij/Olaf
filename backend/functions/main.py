import json
import time
from executor import Executor
from agents.coder_agent import CoderAgent
from agents.master_agent import MasterAgent
import flask

from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
import concurrent.futures
from history import History
from functions_framework import http
from agent_utils import chat_completion, stream
initialize_app()

# import other modules' Cloud Functions
from e2b_functions import (
    request_sandbox,
    execute_on_sandbox,
    upload_to_sandbox,
    download_from_sandbox,
    sandbox_status,
    close_sandbox,
    firebase_storage_to_sandbox
)

from file_storage_functions import (
    handle_user_file_upload,
    request_user_create_folder,
    request_user_delete_path
)

from sessions_functions import (
    update_session,
    delete_session,
    get_sessions,
    delete_all_sessions,
    create_session
)

# This needs to be cleaned up
E2B_API_KEY = "REMOVED"

# --- EXAMPLE ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get", "post"]))
def on_request_example(req: Request) -> Response:
    return Response("Hello world this is me!!")

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def master_agent_interaction(req: Request) -> Response:
    try:
        history = req.json.get("history")
        
        if not history:
            return flask.Response(json.dumps({"error": "'history' is required"}), status=400)
        
        history = History(history)
        master_agent = MasterAgent(history)

        return flask.Response(flask.stream_with_context(stream(master_agent)), mimetype="text/event-stream")
    except Exception as e:
        print(f"Error in generate_plan: {str(e)}")
        return flask.Response(json.dumps({"error": str(e)}), status=500)


# --- CoderAgent Functions ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_code(req: Request) -> flask.Response:
    try:
        history = req.json.get("history")
        language = req.json.get("language", "Python")
        
        if not history:
            return flask.Response(json.dumps({"error": "'history' is required"}), status=400)
        
        history = History(history)
        coder_agent = CoderAgent(language=language, history=history)
    
        return flask.Response(flask.stream_with_context(stream(coder_agent)), mimetype="text/event-stream")
    except Exception as e:
        print(f"Error in generate_code: {str(e)}")
        return flask.Response(json.dumps({"error": str(e)}), status=500)

# --- LLM Agent Function ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def ask_agent(req: Request) -> Response:
    # a request will be conversation history
    system = """You are a highly knowledgeable and experienced expert in the field of bioinformatics. Your primary role is to provide clear, concise answers and general context on topics related to bioinformatics. However, your responses should be limited to offering explanations, insights, and overviews. You are not to provide detailed plans, step-by-step instructions, or complete code. Your expertise serves as a foundation for someone else to build upon, and they will use your responses to generate specific plans and code.

    Guidelines:

        1.	Scope of Responses:
        •	Focus on providing background information, definitions, and explanations.
        •	Offer high-level guidance and theoretical insights.
        •	Clarify complex concepts and discuss relevant methodologies in general terms.
        2.	Limitations:
        •	Avoid giving specific, actionable plans or instructions.
        •	Do not provide or suggest specific code or algorithms.
        •	Refrain from diving into implementation details or technical procedures.
        3.	Purpose:
        •	Your responses will serve as a resource for another individual who will synthesize your insights into actionable plans and code.
        •	Aim to educate and inform, laying the groundwork for further development by others.

    Examples:

        •	Good Response: “In bioinformatics, sequence alignment is a fundamental method used to identify regions of similarity between DNA, RNA, or protein sequences. It’s commonly used in identifying evolutionary relationships and functional similarities between sequences. Tools like BLAST are frequently used for this purpose.”
        •	Inappropriate Response: “To align sequences, you can use the following Python code snippet with the Biopython library: from Bio import pairwise2…”

    This prompt clarifies your role, the scope of your responses, and the limitations of what you should provide, ensuring that you focus on high-level expertise without delving into detailed implementation.
    """
    history = req.json.get("history")
    history = History(history)
    history.log("system", system)
    response = chat_completion(history)
    print(response)
    response_data = {
        "message": response,
    }
    
    return Response(json.dumps(response_data), status=200)
