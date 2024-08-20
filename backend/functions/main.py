import json
from executor import Executor
from agents.coder_agent import CoderAgent
from agents.tester_agent import TesterAgent
from agents.master_agent import MasterAgent

from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

import concurrent.futures
from functions_framework import http

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
        
        master_agent = MasterAgent()
        plan = master_agent.plan(history)
        
        response_data = {
            "plan": plan
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)

@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_plan(req: Request) -> Response:
    try:
        history = req.json.get("history")
        previous_plan = req.json.get("previous_plan")
        
        if not history or not previous_plan:
            return Response(json.dumps({"error": "Both 'history' and 'previous_plan' are required"}), status=400)
        
        master_agent = MasterAgent()
        plan = master_agent.re_plan(history, previous_plan)
        
        response_data = {
            "plan": plan
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- CoderAgent Functions ---
@http
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_code(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        language = req.json.get("language", "Python")
        
        if not plan:
            return Response(json.dumps({"error": "'plan' is required"}), status=400)
        
        coder_agent = CoderAgent(language=language)
        
        generated_code = coder_agent.generate(plan)
        
        response_data = {
            "code": generated_code
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
            "code": regenerated_code
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
            "tests": generated_tests
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
            "tests": regenerated_tests
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
    history = req.json.get("history")
    agent = MasterAgent()
    response = agent.chat_completion(history)
    parsed_response = agent.parse_output(response)

    return Response(parsed_response)
