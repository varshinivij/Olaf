import json

from executor import Executor
from agents.coder_agent import CoderAgent
from agents.tester_agent import TesterAgent
from agents.master_agent import MasterAgent

from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from e2b import Sandbox
from e2b_code_interpreter import CodeInterpreter
import concurrent.futures

# initialize cloud functions defined in this file  (organizing like this
# may be better). maybe even put helper modules in a separate /utils folder.
from file_storage_functions import (
    handle_user_file_upload,
    request_user_create_folder,
    request_user_delete_path,
)


initialize_app()

# Global dictionaries to store agent instances
agent_store = {
    "coder_agent": None,
    "tester_agent": None
}

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


# --- MasterAgent Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def plan_query(req: Request) -> Response:
    query = req.json.get("query")
    master_agent = MasterAgent()
    plan = master_agent.plan(query)
    
    response_data = {
        "plan": plan
    }
    
    return Response(json.dumps(response_data), status=200)


# --- CoderAgent Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_code(req: Request) -> Response:
    plan = req.json.get("plan")
    language = req.json.get("language", "Python")
    coder_agent = CoderAgent(language=language)
    agent_store["coder_agent"] = coder_agent
    
    generated_code = coder_agent.generate(plan)
    
    response_data = {
        "code": generated_code
    }
    
    return Response(json.dumps(response_data), status=200)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_code(req: Request) -> Response:
    plan = req.json.get("plan")
    code_result = req.json.get("code_result")
    test_result = req.json.get("test_result")
    language = req.json.get("language", "Python")
    coder_agent = agent_store.get("coder_agent")
    
    regenerated_code = coder_agent.regenerate(plan, code_result, test_result)
    
    response_data = {
        "code": regenerated_code
    }
    
    return Response(json.dumps(response_data), status=200)


# --- TesterAgent Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_tests(req: Request) -> Response:
    plan = req.json.get("plan")
    language = req.json.get("language", "Python")
    tester_agent = TesterAgent(language=language)
    agent_store["tester_agent"] = tester_agent
    
    generated_tests = tester_agent.generate(plan)
    
    response_data = {
        "tests": generated_tests
    }
    
    return Response(json.dumps(response_data), status=200)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_tests(req: Request) -> Response:
    plan = req.json.get("plan")
    code_result = req.json.get("code_result")
    test_result = req.json.get("test_result")
    language = req.json.get("language", "Python")
    tester_agent = agent_store.get("tester_agent")
    
    regenerated_tests = tester_agent.regenerate(plan, code_result, test_result)
    
    response_data = {
        "tests": regenerated_tests
    }
    
    return Response(json.dumps(response_data), status=200)


# --- Workflow Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def code_generation_workflow(req: Request) -> Response:
    planner_agent = MasterAgent()
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