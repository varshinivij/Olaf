import json
from executor import Executor
from agents.coder_agent import CoderAgent
from agents.master_agent import MasterAgent

from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

import concurrent.futures
from history import History
from functions_framework import http
from agent_utils import chat_completion

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

# This needs to be cleaned up
E2B_API_KEY = "REMOVED"

# --- EXAMPLE ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get", "post"]))
def on_request_example(req: Request) -> Response:
    return Response("Hello world this is me!!")


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_plan(req: Request) -> Response:
    try:
        history_data = req.json.get("history")
        
        if not history_data:
            return Response(json.dumps({"error": "'history' is required"}), status=400)
        
        history = History(history_data)
        master_agent = MasterAgent(history)
        plan = master_agent.process_query(history)
        
        response_data = {
            "message": plan
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- CoderAgent Functions ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_code(req: Request) -> Response:
    try:
        history = req.json.get("history")
        history = History(history)
        language = req.json.get("language", "Python")
        
        if not history:
            return Response(json.dumps({"error": "'history' is required"}), status=400)
        
        coder_agent = CoderAgent(language=language)
        
        generated_code = coder_agent.generate(history)

        response_data = {
            "message": generated_code
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)




# http://127.0.0.1:4000/twocube-web//us-central1/?query=I%20have%20uploaded%2016%20files,%204%20files%20per%20cell-type.%20For%20each%20cell%20type%20there%20are%20three%20negative%20sequence%20files%20and%20one%20positive%20sequence%20file.%20Build%20a%20convolution%20neural%20network%20based%20model%20to%20classify%20positive%20and%20negative%20DNA%20sequences.%20For%20evaluation%20results,%20plot%20the%20area%20under%20precision%20recall%20curve%20and%20area%20under%20the%20receiver%20operator%20characteristic%20curve.
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
