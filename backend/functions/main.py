import json
from executor import Executor
from agents.coder_agent import CoderAgent
from agents.tester_agent import TesterAgent
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

# --- MasterAgent Functions ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_plan(req: Request) -> Response:
    try:
        history = req.json.get("history")
        if not history:
            return Response(json.dumps({"error": "'history' is required"}), status=400)
        history = History(history)
        master_agent = MasterAgent()
        plan = master_agent.plan(history)
        
        response_data = {
            "message": plan
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_plan(req: Request) -> Response:
    try:
        history = req.json.get("history")
        history = History(history)
        previous_plan = req.json.get("previous_plan")
        
        if not history or not previous_plan:
            return Response(json.dumps({"error": "Both 'history' and 'previous_plan' are required"}), status=400)
        
        master_agent = MasterAgent()
        plan = master_agent.re_plan(history, previous_plan)
        
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

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_code(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        code_result = req.json.get("code_result")
        test_result = req.json.get("test_result")
        language = req.json.get("language", "Python")
        
        if not plan or not code_result or not test_result:
            return Response(json.dumps({"error": "All of 'plan', 'code_result', and 'test_result' are required"}), status=400)
        
        coder_agent = CoderAgent(language=language)
        
        regenerated_code = coder_agent.regenerate(plan, code_result, test_result)
        
        response_data = {
            "message": regenerated_code
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- TesterAgent Functions ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_tests(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        language = req.json.get("language", "Python")
        
        if not plan:
            return Response(json.dumps({"error": "'plan' is required"}), status=400)
        
        tester_agent = TesterAgent(language=language)
        
        generated_tests = tester_agent.generate(plan)
        
        response_data = {
            "message": generated_tests
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_tests(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        code_result = req.json.get("code_result")
        test_result = req.json.get("test_result")
        language = req.json.get("language", "Python")
        
        if not plan or not code_result or not test_result:
            return Response(json.dumps({"error": "All of 'plan', 'code_result', and 'test_result' are required"}), status=400)
        
        tester_agent = TesterAgent(language=language)
        
        regenerated_tests = tester_agent.regenerate(plan, code_result, test_result)
        
        response_data = {
            "message": regenerated_tests
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- Workflow Function ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def code_generation_workflow(req: Request) -> Response:
    try:
        history = req.json.get("history")
        if not history:
            return Response(json.dumps({"error": "'history' is required"}), status=400)
        history = History(history)
        planner_agent = MasterAgent()
        coder_agent = CoderAgent(language='Python')
        tester_agent = TesterAgent(language='Python')
        executor_agent = Executor(api_key=E2B_API_KEY)
        
        requirements = planner_agent.plan(history)
        
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_code = executor.submit(coder_agent.generate, requirements)
            future_tests = executor.submit(tester_agent.generate, requirements)
            
            generated_code = future_code.result()
            generated_tests = future_tests.result()

        count = 0
        code_result = "Failed"
        tests_result = "Failed"
        
        while (code_result != "Passed" or tests_result != "Passed") and count < 5:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future_code_result = executor.submit(executor_agent.execute_code, {"code": generated_code})
                future_tests_result = executor.submit(executor_agent.execute_tests, {"code": generated_code, "tests": generated_tests})
                
                code_result = future_code_result.result().json.get("output", "Failed")
                tests_result = future_tests_result.result().json.get("output", "Failed")

            if code_result != "Passed" or tests_result != "Passed":
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future_code = executor.submit(coder_agent.regenerate, requirements, code_result, tests_result)
                    future_tests = executor.submit(tester_agent.regenerate, requirements, code_result, tests_result)
                    
                    generated_code = future_code.result()
                    generated_tests = future_tests.result()

            count += 1
            if count == 4:
                return Response(json.dumps({"error": "Iteration limit exceeded :("}), status=500)
        
        response_data = {
            "code": generated_code,
            "tests": generated_tests,
            "code_result": code_result,
            "tests_result": tests_result
        }
        
        return Response(json.dumps(response_data), status=200)

    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)
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
