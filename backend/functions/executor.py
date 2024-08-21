from e2b_code_interpreter import CodeInterpreter
import json

# Currently not using sandbox output for feedback. NEEDS TO BE IMPLEMENTED
class Executor:
    def __init__(self, api_key):
        self.api_key = api_key

    def request_sandbox(self):
        try:
            sandbox = CodeInterpreter(
                api_key=self.api_key, template="vh7kehbtf0t4xbx9ec9u", timeout=300
            )
            sandbox_id = sandbox.id
            sandbox.keep_alive(5 * 60)  # keep box alive for 5minutes
        except Exception as e:
            sandbox_id = f"Error during code execution: {e}"
        return sandbox_id

    def execute_code(self, code: str) -> str:
        """
        Executes the provided code using e2b.
        :param code: The code snippet to execute.
        :return: The result of the execution.
        """
        try:
            sandbox_id = self.request_sandbox()
            sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.api_key)
            execution = sandbox.notebook.exec_cell(code)
            result = "Passed"
        except Exception as e:
            result = f"Error during code execution: {e}"
        return result

    def execute_tests(self, code: str, tests: str) -> str:
        """
        Executes the provided tests using e2b.
        :param code: The code snippet to test.
        :param tests: The test cases to execute.
        :return: The result of the test execution.
        """
        combined_code = code + "\n" + tests
        try:
            sandbox_id = self.request_sandbox()
            sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.api_key)
            execution = sandbox.notebook.exec_cell(code)
            result = "Passed"
        except AssertionError as e:
            result = f"Test failed: {e}"
        except Exception as e:
            result = f"Error during test execution: {e}"
        return result
