from functions.agents.coder_agent import CoderAgent
from functions.agents.tester_agent import TesterAgent
from functions.agents.planner_agent import PlannerAgent
from executor import Executor
from firebase_admin import initialize_app

from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from e2b import Sandbox
from e2b_code_interpreter import CodeInterpreter
import concurrent.futures

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


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def code_generation_workflow(req: Request) -> Response:
    planner_agent = PlannerAgent()
    coder_agent = CoderAgent(language='Python')
    tester_agent = TesterAgent(language='Python')
    executor_agent = Executor(api_key=E2B_API_KEY)
    
    query = req.json.get("query")

    requirements = planner_agent.plan(query)
    
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
            return Response("Iteration limit exceeded :(", status=500)
    
    response_data = {
        "code": generated_code,
        "tests": generated_tests,
        "code_result": code_result,
        "tests_result": tests_result
    }
    
    return Response(json.dumps(response_data), status=200)
