from prompts import few_shot_examples
from agent_utils import *
from history import History

system_prompt = """
You are `CoderAgent`, a programming assistant specialized in generating and refining code. Your task is to:
- Generate code in the specified language based on provided requirements.
- Revise and improve code using feedback from code results and test outcomes.
- Exclude test cases from the final code output during refinement.
"""

system = {
        "role": "system",
        "content": system_prompt
      }


class CoderAgent:

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
        prompt = f"{few_shot_examples}\n\n\nGenerate {self.language} program that meets the specified requirements:\n{requirements}"
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
        {few_shot_examples}\n\n\n
        Based on the results provided below, please revise and improve the program to ensure it functions correctly. 
        Generate a new {self.language} program that meets the specified requirements:\n{requirements}.
        Do not generate the tests itself in the output, only the {self.language} program.

        code_result: {code_result}
        test_result: {test_result}
        """
        self.history.log("user", prompt)
        response = chat_completion(self.history)
        self.history.log("assistant", response)
        return extract_python_code(response)