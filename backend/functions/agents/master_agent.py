import openai
from agent_utils import chat_completion_api
from agent_utils import chat_completion_plan
from agent_utils import chat_completion_function
import json
from firebase_functions.https_fn import Request, Response, on_request

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
        - You are a knowledgeable assistant. Please provide a clear and concise response 
          Your response should be easy to understand and directly address the query.


    2. Basic Code Writing: -> In this case please use write_basic_code function call please 
        - Write simple code snippets when necessary for straightforward tasks.
        - when simple code tasks are mentioned this is done
        - please create like basic code when needed


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
        
        
        
        f"You are an expert in project planning, especially in the domains of bioinformatics, machine learning, and software development. "
        f"The user has provided a complex query that requires a step-by-step plan to execute. "
        "Please create a comprehensive sequential plan that outlines each necessary step clearly and logically. "
        "Ensure that the plan is detailed enough for implementation, and include any dependencies or prerequisites required for each step."


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
    def __init__(self,history):
        self.system_prompt = system_prompt
        self.history = history
        self.functions = [
            {
                "type": "function",
                "function": {
                    "name": "handle_simple_interaction",
                    "description": "Handles simple user interactions. This is for like human interactions or something like that.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                        },
                        "required": ["history"],
                        "additionalProperties": False
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "write_basic_code",
                    "description": "This function writes simple python code queries where not a lot of computations or code is required use this when going for very small functions . Use this function for creating code when very small piece of code is asked like simple functions code or simple class or anything.",
                    "parameters": {
                        "type": "object",
                        "properties": {
                        },
                        "required": ["history"],
                        "additionalProperties": False
                    }
                }
            },
            {
                "type": "function",
                "function": {
                    "name": "create_sequential_plan",
                    "description": "Creates a detailed, step-by-step plan for a complex task or larger code-related query. Use this function for planning out intricate workflows or comprehensive code architectures. This function should also be used if asked for regenrate or modify the plan okay. ",
                    "parameters": {
                        "type": "object",
                        "properties": {
                        },
                        "required": ["history"],
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

    def generate(self):
        try:
            self.history.remove_system_messages()
            self.history.upsert(self.system_prompt)
            # Call the chat completion API
            response = chat_completion_function(self.history, tools=self.functions)
            
            # Check if response is a dictionary (API response)
            if isinstance(response, openai.types.chat.chat_completion.ChatCompletion):
                # Extract the message from the API response
                message = response.choices[0].message

                # Check if there are tool calls
                if message.tool_calls:
                    # Handle tool calls
                    tool_call = message.tool_calls[0]
                    function_name = tool_call.function.name
                    function_arguments = tool_call.function.arguments
                    print(function_name)
                    # Ensure function name is present
                    if not function_name:
                        return {"error": "Function name is missing in the tool call."}

                    # Convert the arguments from JSON to a Python dictionary if they exist
                    args_dict = {}
                    if function_arguments:
                        try:
                            args_dict = json.loads(function_arguments)
                        except json.JSONDecodeError:
                            return {"error": "Function arguments could not be parsed as JSON."}

                    # Check if the function exists in the function map
                    if function_name in self.function_map:
                        print("function name : ",function_name)
                        try:
                            # Call the appropriate function with or without arguments
                            if args_dict:
                                result = self.function_map[function_name](**args_dict)
                            else:
                                result = self.function_map[function_name]()
                            return result
                        except Exception as e:
                            return {"error": f"Error calling function {function_name}: {str(e)}"}
                    else:
                        return {"error": f"Function {function_name} not found in function map."}
                else:
                    # Return the content if no tool calls
                    return response.get('content', '')
            else:
                # If response is not a dictionary, assume it's the direct content
                return response
        except json.JSONDecodeError:
            return {"error": "Invalid response from API (JSON decode error)."}
        except Exception as e:
            return {"error": f"An unexpected error occurred: {str(e)}"}

    def handle_simple_interaction(self):
        interaction_response = chat_completion_api(self.history, system_prompt)
        result = ""
        for chunk in interaction_response:
            try:
                content = chunk['choices'][0]['delta']['content']
                if content:
                    result += content
                yield chunk
            except:
                continue
        yield "Response: text"   
        self.history.log("assistant", result)

    def write_basic_code(self):
        result = ""
        for chunk in chat_completion_api(self.history, system_prompt):
            try:
                print(chunk)
                content = chunk['choices'][0]['delta']['content']
                if content:
                    result += content
                yield chunk
            except:
                continue
        yield "Response: code"   
        self.history.log("assistant", result)

    def decompose_complicated_task(self):
        decomposition_prompt = (
            f"You are an expert in task decomposition, particularly in bioinformatics and machine learning. "
            f"The user has provided a complex query that needs to be broken down into smaller, manageable subtasks. "
            f"Here is the query:\n\n"
            f"Query: '{self.history.most_recent_entry()}'\n\n"
            "Please identify the key components of this task and provide a step-by-step breakdown into smaller subtasks."
        )
        self.history.remove_system_messages()
        self.history.log("user", decomposition_prompt)

        result = ""
        for chunk in chat_completion_api(self.history, system_prompt):
            try:
                content = chunk['choices'][0]['delta']['content']
                if content:
                    result += content
                yield chunk
            except:
                continue
        yield "Response: task"   
        self.history.log("assistant", result)

    def create_sequential_plan(self):
        self.history.remove_system_messages()
        result = ""
        for chunk in chat_completion_plan(self.history, system_prompt):
            try:
                content = chunk['choices'][0]['delta']['content']
                if content:
                    result += content
                yield chunk
            except:
                continue
        yield "Response: plan"    
        self.history.log("assistant", result)
    
        

