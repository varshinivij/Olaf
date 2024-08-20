# from history import History
# from agent_utils import chat_completion
# from agent_utils import chat_completion_api
import json
from openai import OpenAI
from pydantic import BaseModel
import requests
import re


class History:
    def __init__(self):
        self.entries = []

    def add_entry(self, role, content):
        self.entries.append({"role": role, "content": content})

    def get_history(self):
        return self.entries

    def clear_history(self):
        self.entries = []

    def __repr__(self):
        return "\n".join([f"{entry['role'].capitalize()}: {entry['content']}" for entry in self.entries])


system_prompt = """
    You are MasterAgent, an AI assistant specialized in bioinformatics. Your responsibilities include:

    1. User Interaction: whenever its just simple user interaction please use the function -> handle_simple_interaction . when simply user is trying to have a conversation . 
        - Engage in regular conversations with users.
        - Collect detailed descriptions of data and tasks.
        - Ask follow-up questions wherever required to ensure comprehensive understanding.
        - Like if its normal greeting or queries where everything can be answered in one go there is no requirement of followup but if query require follow up , you can ask follow up questions.
        - If its only a greeting there is no need to classify just greet them back please 
        - if its not realted to bioinformatics or software or something you can make your punchline . 
        - if anything else is asked you can process the query and respond.

    2. Basic Code Writing:
        - Write simple code snippets when necessary for straightforward tasks.
        - when simple code tasks are mentioned this is done


        Examples:
        - "Write a Python function to calculate the sum of two numbers." -> (Write the function directly)
        - "Show me how to reverse a string in Python." -> (Provide the code for string reversal)
        - "How do I read a CSV file in Python?" -> (Provide a simple code example using `pandas`)
        - "Create a function to multiply two matrices in Python." -> (Write the matrix multiplication code)
        - "How do I generate a random number in Python?" -> (Provide the code using the `random` module)
        - "Write a function to convert Fahrenheit to Celsius." -> (Write the conversion function)
        - "Can you create a list comprehension to filter even numbers?" -> (Provide the list comprehension code)
        - "How do I write to a text file in Python?" -> (Provide the code to write to a file)
        - "Write a Python function to check if a number is prime." -> (Provide the prime-checking code)
        - "Create a simple `for` loop in Python." -> (Write the `for` loop example)


    3. Sequential Planning: -> if larger code queries or any sort ofneural network or anything is asked you have to use the create_sequential_plan function please . Anything that is not related to basic code execution and simple convo you
    have to do sequential planning as the function call please 
        - This is like when the planning is required for complex tasks ok  
        - Plan the sequence of actions to achieve a goal.
        - Generate detailed step-by-step plans for complex projects.
        - Display these plans to the user after all information is collected.
        - Till the user is not satisfied with the plan please ask follow up questions like do we need to change 
        - This is also when a larger query is given like a larger code query is given so we need to make a plan before actually coding it.


        Examples:
        - "Develop a neural network for protein classification." -> (Create a plan outlining data preparation, model design, training, validation, etc.)
        - "Build a pipeline for whole-genome sequencing analysis." -> (Plan out data acquisition, alignment, variant calling, annotation, etc.)
        - "Design a cloud-based bioinformatics platform." -> (Create a step-by-step plan for infrastructure setup, security, scalability, etc.)
        - "Automate the entire RNA-Seq workflow." -> (Plan out read alignment, quantification, differential expression, visualization, etc.)
        - "Implement a machine learning model for predicting drug resistance." -> (Plan data collection, feature engineering, model selection, evaluation, etc.)
        - "Set up a database for large-scale genomic data." -> (Create a plan for schema design, data storage, indexing, query optimization, etc.)
        - "Design a tool for single-cell RNA-Seq analysis." -> (Plan out preprocessing, clustering, differential expression, visualization, etc.)
        - "Create a bioinformatics pipeline for metagenomics." -> (Plan out sample processing, sequencing, taxonomic classification, functional analysis, etc.)
        - "Develop a visualization tool for proteomics data." -> (Create a plan for data integration, visual design, interactive features, etc.)
        - "Implement a CRISPR screening analysis workflow." -> (Plan out guide RNA design, screening, analysis, validation, etc.)
        - I have provided files please build a neural network or build a code to analyse the files and give analysis . 
        - if the user asks to modify the plan you should again run this with the previous plan and then the user query okay please 
        


    5. Iterative Optimization:
        - Refine and optimize plans based on user feedback.

    6. Coder Agent Interaction:
        - Coordinate with the CoderAgent for tasks requiring extensive code generation or technical implementation.


    points : 

    1. Make plans only if its a big tasks or complex tasks not for smaller interactions and basic code.
    2. Use the supplied tools to assist the user.I will provide you with some tools which you need to use to assist the user
    3. MOST IMPORTANT : ANY LARGER QUERIES OR LARGER CODE GENERATION FIRST YOU NEED TO GENERATE A PLAN USING THE CREATE SEQUENTIAL PLAN FUNCTIONS
  
    """


system = {
    "role": "system",
    "content": system_prompt
}


class MasterAgent:
    def __init__(self, system_prompt):
        self.client = OpenAI(
            api_key="REMOVED")
        self.system_prompt = system_prompt
        self.history = History()
        self.functions = [
            {
                "type": "function",
                "function": {
                    "name": "handle_simple_interaction",
                    "description": "Handles simple user interactions",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "query": {
                                "type": "string",
                                "description": "The user's query"
                            }
                        },
                        "required": ["query"],
                        "additionalProperties": False
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "write_basic_code",
                    "description": "This function writes simple code queries where not a lot of computations or code is required",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "query": {
                                "type": "string",
                                "description": "The user's query"
                            }
                        },
                        "required": ["query"],
                        "additionalProperties": False
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "create_sequential_plan",
                    "description": "Creates a detailed, step-by-step plan for a complex task or larger code-related query. Use this function for planning out intricate workflows or comprehensive code architectures. This function should also be used if asked for regenrate or modify the plan okay . use the history for previou plan and then use then use user query to generate new plan ",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "query": {
                                "type": "string",
                                "description": "The user's query"
                            }
                        },
                        "required": ["query"],
                        "additionalProperties": False
                    }
                }
            },
        ]

        self.function_map = {
            "handle_simple_interaction": self.handle_simple_interaction,
            "write_basic_code": self.write_basic_code,
            "create_sequential_plan": self.create_sequential_plan,
        }

    def chat_completion_api(self, query, tools=None):
        print(query)
        api_key = "REMOVED"
        client = OpenAI(api_key=api_key)
        completion = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                # Include the history in the conversation
                {"role": "system", "content": self.system_prompt},
                *self.history.get_history(),
                {"role": "user", "content": query}
            ],
            tools=tools
        )
        return completion.choices[0].message

    def process_query(self, query: str):
        try:
            response = self.chat_completion_api(query, tools=self.functions)
            if response.tool_calls:
                # Step 2: Extract the function name and arguments
                # Assuming there's only one tool call
                tool_call = response.tool_calls[0]
                function_name = tool_call.function.name
                function_arguments = tool_call.function.arguments

                # Step 3: Convert the arguments from JSON to a Python dictionary
                args_dict = json.loads(function_arguments)


                # Step 4: Call the appropriate function
                if function_name in self.function_map:
                    result = self.function_map[function_name](**args_dict)
                    return result.content
                else:
                    return {"error": f"Function {function_name} not found."}
            else:
                return response.content
        except json.JSONDecodeError:
            return {"error": "Invalid response from API"}

    def handle_simple_interaction(self, query):
        interaction_prompt = (
            f"You are a knowledgeable assistant. Please provide a clear and concise response "
            f"to the following query:\n\n"
            f"Query: '{query}'\n\n"
            "Your response should be easy to understand and directly address the query."
        )
        self.history.add_entry("user", interaction_prompt)
        interaction_response = self.chat_completion_api(interaction_prompt)
        return interaction_response

    def write_basic_code(self, query):
        code_prompt = (
            f"You are an expert software developer. Based on the user's query, please generate a simple, functional code snippet.\n\n"
            f"Query: '{query}'\n\n"
            "Please write the code below:\n"
        )
        self.history.add_entry("user", code_prompt)
        code_response = self.chat_completion_api(code_prompt)
        self.history.add_entry("assistant", code_response)
        return code_response.strip()

    def decompose_complicated_task(self, query):
        decomposition_prompt = (
            f"You are an expert in task decomposition, particularly in bioinformatics and machine learning. "
            f"The user has provided a complex query that needs to be broken down into smaller, manageable subtasks. "
            f"Here is the query:\n\n"
            f"Query: '{query}'\n\n"
            "Please identify the key components of this task and provide a step-by-step breakdown into smaller subtasks."
        )
        self.history.add_entry("user", decomposition_prompt)
        decomposition_response = self.chat_completion_api(decomposition_prompt)
        return decomposition_response.strip()

    def create_sequential_plan(self, query):
        plan_creation_prompt = (
            f"You are an expert in project planning, especially in the domains of bioinformatics, machine learning, and software development. "
            f"The user has provided a complex query that requires a step-by-step plan to execute. "
            f"Here is the query:\n\n"
            f"Query: '{query}'\n\n"
            "Please create a comprehensive sequential plan that outlines each necessary step clearly and logically. "
            "Ensure that the plan is detailed enough for implementation, and include any dependencies or prerequisites required for each step."
        )
        self.history.add_entry("user", plan_creation_prompt)
        # self.chat_completion_api(plan_creation_prompt)
        plan_response = self.client.chat.completions.create(
            model="gpt-4o-2024-08-06",
            messages=[
                {"role": "system", "content": self.system_prompt},
                *self.history.get_history(),
                {"role": "user", "content": query}
            ],
            response_format={
                "type": "json_schema",
                "json_schema": {
                    "name": "plan_response",
                    "schema": {
                        "type": "object",
                        "properties": {
                            "overview": {
                                "type": "string",
                                "description": "A brief overview of the entire plan."
                            },
                            "steps": {
                                "type": "array",
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "step_number": {"type": "integer"},
                                        "title": {"type": "string", "description": "A short title or name for the step."},
                                        "description": {"type": "string", "description": "A detailed explanation of what this step involves."},
                                        "dependencies": {
                                            "type": "array",
                                            "items": {"type": "string"},
                                            "description": "Any prerequisites or dependencies required for this step."
                                        },
                                        "expected_output": {"type": "string", "description": "The expected result or output of this step."}
                                    },
                                    "required": ["step_number", "title", "description", "dependencies", "expected_output"],
                                    "additionalProperties": False
                                }
                            },
                            "final_check": {
                                "type": "string",
                                "description": "A final review or checklist after all steps are completed."
                            }
                        },
                        "required": ["overview", "steps", "final_check"],
                        "additionalProperties": False
                    },
                    "strict": True
                }
            }
        )
        self.history.add_entry(
            "assistant", plan_response.choices[0].message.content)
        return plan_response.choices[0].message


def main():
    agent = MasterAgent(system_prompt)
    print("Welcome! You can start chatting with the Master Agent.")
    print("Type 'exit' to end the conversation.")

    while True:
        user_input = input("You: ")
        if user_input.lower() == "exit":
            print("Ending conversation. Goodbye!")
            break

        # Process the user query with the Master Agent
        response = agent.process_query(user_input)
        print(response)

        # Display the response
        print(f"Agent: {response}")


if __name__ == "__main__":
    main()
