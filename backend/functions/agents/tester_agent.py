from prompts import few_shot_examples
from functions.agent_utils import *
from history import History

system_prompt = f"""
You are `TesterAgent`, an expert in generating and refining test cases for programming code. Your task is to:
- Generate test cases based on the provided requirements, ensuring all necessary test scenarios are covered.
- Regenerate and improve test cases using feedback from code results and test outcomes.
- Exclude the program code itself from the output, focusing solely on producing comprehensive and accurate test cases.
"""

system = {
        "role": "system",
        "content": system_prompt
      }


class TesterAgent:

    def __init__(self, language: str):
        """
        Initializes the CoderAgent with language and LLM instance.
        :param language: The programming language for code generation.
        """
        self.language = language
        self.history = History(system)

    def generate(self, requirements: str) -> str:
        """
        Generates code based on the provided requirements using the LLM.
        :param requirements: A dictionary containing the requirements for the code.
        :return: A string of the generated code.
        """
        prompt = f"Generate {self.language} tests for a {self.language} program that implements the following requirements:\n{requirements}. Do not generate the program itself in the output, only the {self.language} tests."
        self.history.log("user", prompt)
        response = chat_completion(self.history)
        self.history.log("assistant", response)
        return extract_python_code(response)
    
    def regenerate(self, requirements: str, code_result: str, test_result: str) -> str:
        """
        Generates test cases for the provided code.
        :param code: The code snippet to generate tests for.
        :return: A string of generated test cases.
        """
        prompt = f"""
        Based on the results below, please revise and improve the test cases for the program.
        Generate again only the {self.language} tests for a {self.language} program that implements the following requirements:\n{requirements}.
        Do not generate the program itself in the output, only the {self.language} tests.

        code_result: {code_result}
        test_result: {test_result}
        """
        self.history.log("user", prompt)
        response = chat_completion(self.history)
        self.history.log("assistant", response)
        return extract_python_code(response)