import json
import time
from executor import Executor
from agents.coder_agent import CoderAgent
from agents.master_agent import MasterAgent
from agents.codemaster_agent import CodeMasterAgent
import flask

from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions
import concurrent.futures
from history import History
from functions_framework import http
from agent_utils import chat_completion, stream
from firebase_admin import firestore
initialize_app()
import sessions_functions
from router import Router
from pipe import Pipe

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
    add_message_to_session,
    delete_session,
    get_sessions,
    delete_all_sessions,
    rename_session
)

# This needs to be cleaned up
E2B_API_KEY = "REMOVED"
db = firestore.client()


# --- EXAMPLE ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["get", "post"]))
def on_request_example(req: Request) -> Response:
    return Response("Hello world this is me!!")


def master_route_function(session):
    history = session
    if not history:
            return None
    history = History(history)
    master_agent = CodeMasterAgent("Python", history) #using CodeMasterAgent
    return master_agent
    

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def master_agent_interaction(req: Request) -> Response:
    try:
        session_id = req.args.get("session_id")
        message = req.json.get("message")
        user_id = req.args.get("user_id")
        project_id = req.args.get("project_id")
        
        if not message:
            return Response(json.dumps({"error": "'message' is required"}), status=400)
        if not user_id:
            return Response(json.dumps({"error": "'user_id' is required"}), status=400)
        if not project_id:
            return Response(json.dumps({"error": "'project_id' is required"}), status=400)

        if session_id:
            # Append the new user message to the existing session history
            session_data = sessions_functions.add_message_to_session(session_id, message, role="user")
        else:
            session_id, session_data = sessions_functions.create_new_session(user_id, project_id, message)

        # Route the session history for response generation
        router = Router()
        router.add_route("master", master_route_function)
        updated_session = router.route_session("master", session_data["history"])

        def generate_stream():
            yield session_id.encode('utf-8')
            full_response = ""
            response_type = ""
            for line in stream(updated_session):
                if response_type == "":
                    response_type = line.decode('utf-8')
                else:
                    full_response += line.decode('utf-8')
                yield line
            
            # Update the session history with the complete assistant response
            add_message_to_session(session_id, full_response, role="assistant")

        if updated_session is None:
            return sesponse(json.dumps({"error": "'history' is required"}), status=400)

        return Response(flask.stream_with_context(generate_stream()), headers={"session_id": session_id}, mimetype="text/event-stream")
    
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


@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def name_maker(req: Request) -> Response:
    # a request will be conversation history
    system = """
    You are given a conversation history.
    Please construct a super short title for the conversation.
    """
    history = req.json.get("history")
    history = History(history)
    history.log("system", system)
    response = chat_completion(history)
    # Initialize an empty string to accumulate the response text
    response_str = ""
    for chunk in response:
        response_str += chunk.get("choices", [{}])[0].get("delta", {}).get("content", "")
    response_data = {
        "message": response_str,
    }
    # Return the response as a JSON response
    return Response(json.dumps(response_data), status=200)

