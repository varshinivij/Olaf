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

# Initialize Firebase app
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


# --- E2B Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def request_sandbox(req: Request) -> Response:
    try:
        sandbox = CodeInterpreter(
            api_key=E2B_API_KEY, template="vh7kehbtf0t4xbx9ec9u", timeout=300
        )
        sandbox_id = sandbox.id
        sandbox.keep_alive(5 * 60)  # keep box alive for 5 minutes
        result = {"sandboxId": sandbox_id}
        json_result = json.dumps(result)

        return Response(json_result, status=200)

    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def execute_on_sandbox(req: Request) -> Response:
    try:
        sandbox_id = req.json.get("sandboxId")
        code = req.json.get("code")

        if not sandbox_id or not code:
            return Response(json.dumps({"error": "Both 'sandboxId' and 'code' are required"}), status=400)

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=E2B_API_KEY)
        execution = sandbox.notebook.exec_cell(code)

        result = {"output": execution.logs.stdout}
        json_result = json.dumps(result)

        return Response(json_result, status=200)

    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_plan(req: Request) -> Response:
    try:
        query = req.json.get("query")
        if not query:
            return Response(json.dumps({"error": "'query' is required"}), status=400)
        
        master_agent = MasterAgent()
        # # plan = master_agent.plan(query)
        # decision = master_agent.make_decision(query)
        # response = ""
        # if(decision=="1"):
        #     response = master_agent.handle_simple_interaction(query)
        # elif(decision=="2"):
        #     response = master_agent.write_basic_code(query)
        # elif(decision=="3"):
        #     response = master_agent.decompose_complicated_task(query)
        # elif(decision=="4"):
        response = master_agent.create_sequential_plan(query)
        
        response_data = {
            "response": response
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_plan(req: Request) -> Response:
    try:
        query = req.json.get("query")
        previous_plan = req.json.get("previous_plan")
        
        if not query or not previous_plan:
            return Response(json.dumps({"error": "Both 'query' and 'previous_plan' are required"}), status=400)
        
        master_agent = MasterAgent()
        plan = master_agent.re_plan(query, previous_plan)
        
        response_data = {
            "plan": plan
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- CoderAgent Functions ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_code(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        language = req.json.get("language", "Python")
        
        if not plan:
            return Response(json.dumps({"error": "'plan' is required"}), status=400)
        
        coder_agent = CoderAgent(language=language)
        agent_store["coder_agent"] = coder_agent
        
        generated_code = coder_agent.generate(plan)
        
        response_data = {
            "code": generated_code
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def regenerate_code(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        code_result = req.json.get("code_result")
        test_result = req.json.get("test_result")
        language = req.json.get("language", "Python")
        
        if not plan or not code_result or not test_result:
            return Response(json.dumps({"error": "All of 'plan', 'code_result', and 'test_result' are required"}), status=400)
        
        coder_agent = agent_store.get("coder_agent")
        
        regenerated_code = coder_agent.regenerate(plan, code_result, test_result)
        
        response_data = {
            "code": regenerated_code
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- TesterAgent Functions ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def generate_tests(req: Request) -> Response:
    try:
        plan = req.json.get("plan")
        language = req.json.get("language", "Python")
        
        if not plan:
            return Response(json.dumps({"error": "'plan' is required"}), status=400)
        
        tester_agent = TesterAgent(language=language)
        agent_store["tester_agent"] = tester_agent
        
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
        
        tester_agent = agent_store.get("tester_agent")
        
        regenerated_tests = tester_agent.regenerate(plan, code_result, test_result)
        
        response_data = {
            "tests": regenerated_tests
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)


# --- Workflow Function ---
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def code_generation_workflow(req: Request) -> Response:
    try:
        query = req.json.get("query")
        if not query:
            return Response(json.dumps({"error": "'query' is required"}), status=400)

        planner_agent = MasterAgent()
        coder_agent = CoderAgent(language='Python')
        tester_agent = TesterAgent(language='Python')
        executor_agent = Executor(api_key=E2B_API_KEY)
        
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
    




@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def simple_convo(req: Request) -> Response:
    try:
        query = req.json.get("query")
        if not query:
            return Response(json.dumps({"error": "'query' is required"}), status=400)
        
        master_agent = MasterAgent()
        response = master_agent.process_query(query)
        
        response_data = {
            "response": response
        }
        
        return Response(json.dumps(response_data), status=200)
    
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)




# http://127.0.0.1:4000/twocube-web//us-central1/?query=I%20have%20uploaded%2016%20files,%204%20files%20per%20cell-type.%20For%20each%20cell%20type%20there%20are%20three%20negative%20sequence%20files%20and%20one%20positive%20sequence%20file.%20Build%20a%20convolution%20neural%20network%20based%20model%20to%20classify%20positive%20and%20negative%20DNA%20sequences.%20For%20evaluation%20results,%20plot%20the%20area%20under%20precision%20recall%20curve%20and%20area%20under%20the%20receiver%20operator%20characteristic%20curve.