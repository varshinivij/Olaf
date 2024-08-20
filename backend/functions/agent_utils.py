import requests
import re
from openai import OpenAI


def chat_completion(history, tools=None):
    api_key = "REMOVED"
    url = "https://api.openai.com/v1/chat/completions"

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }

    payload = {
        "model": "gpt-4o",
        "messages": history.get_history(),
        "temperature": 0.1,
        "tools":  tools
    }

    response = requests.post(url, headers=headers, json=payload)

    if response.status_code == 200:
        result = response.json()
        response_message = result['choices'][0]['message']['content']
        return response_message

    else:
        print(f"Error: {response.status_code}, {response.text}")
        return None


def extract_python_code(text):
    """
    Extract the Python code from the provided text string.

    Args:
    text (str): The input text containing the Python code

    Returns:
    str: The extracted Python code
    """
    code_pattern = re.compile(r'```python(.*?)```', re.DOTALL)
    match = code_pattern.search(text)

    if match:
        return match.group(1).strip()
    else:
        return "No Python code found in the input text."




import json



system_prompt = """
    You are MasterAgent, an AI assistant specialized in bioinformatics. Your responsibilities include:

    1. User Interaction:
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


    4. Sequential Planning:
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


    5. Iterative Optimization:
        - Refine and optimize plans based on user feedback.

    6. Coder Agent Interaction:
        - Coordinate with the CoderAgent for tasks requiring extensive code generation or technical implementation.


    points : 

    1. Make plans only if its a big tasks or complex tasks not for smaller interactions and basic code.
    2. Use the supplied tools to assist the user.I will provide you with some tools which you need to use to assist the user
    3. MOST IMPORTANT : ANY LARGER QUERIES OR LARGER CODE GENERATION FIRST YOU NEED TO GENERATE A PLAN USING THE CREATE SEQUENTIAL PLAN FUNCTIONS
  
    """
tools = [
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
            "description": "Creates a detailed, step-by-step plan for a complex task or larger code-related query. Use this function for planning out intricate workflows or comprehensive code architectures.",
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
            "name": "regenerate_plan",
            "description": "Regenerates and updates an existing plan based on user feedback and new requirements",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "The user's new request or feedback"
                    },
                    "previous_plan": {
                        "type": "string",
                        "description": "The previously generated plan"
                    }
                },
                "required": ["query", "previous_plan"],
                "additionalProperties": False
            }
        }
    }
]
def chat_completion_api(query,system,tools):
    api_key = "REMOVED"
    client = OpenAI(api_key=api_key)
    completion = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": system},
            {
                "role": "user",
                "content": query
            }
        ],
        tools=tools
    )
    #return "Hello"
    # print("message :  " , completion.choices[0].message)
    return completion.choices[0].message

def process_query(query: str, functions, system_prompt):
    try:
        function_call = chat_completion_api(query, system_prompt, functions)
        return function_call
    except json.JSONDecodeError:
        return {"error": "Invalid response from API"}

def handle_simple_interaction(query):
    """
    Handles simple user interactions by generating a concise and clear response 
    to straightforward or non-technical queries.
    """
    interaction_prompt = (
        f"You are a knowledgeable assistant. Please provide a clear and concise response "
        f"to the following query:\n\n"
        f"Query: '{query}'\n\n"
        "Your response should be easy to understand and directly address the query."
    )

    interaction_response = chat_completion_api(interaction_prompt)
    return interaction_response + "hahahahaha"

def write_basic_code(query):
    """
    Handles the creation of simple code for the user's request using the LLM.
    The code generated will be basic and should not require complex logic or architecture.
    """
    code_prompt = (
        f"You are an expert software developer. Based on the user's query, please generate a simple, functional code snippet.\n\n"
        f"Query: '{query}'\n\n"
        "Please write the code below:\n"
    )

    code_response = chat_completion_api(code_prompt,system_prompt,tools=None)
    return code_response

def decompose_complicated_task(query):
    """
    Decomposes a complex query into smaller, manageable subtasks.

    This function identifies the key components of a complicated task and breaks them down
    into a sequence of smaller steps that can be addressed individually.
    """
    decomposition_prompt = (
        f"You are an expert in task decomposition, particularly in bioinformatics and machine learning. "
        f"The user has provided a complex query that needs to be broken down into smaller, manageable subtasks. "
        f"Here is the query:\n\n"
        f"Query: '{query}'\n\n"
        "Please identify the key components of this task and provide a step-by-step breakdown into smaller subtasks."
    )

    decomposition_response = chat_completion_api(decomposition_prompt)
    return decomposition_response.strip() + "hdhhdhdhhdjdjskksk"

def create_sequential_plan(query):
    """
    Creates a detailed, step-by-step sequential plan for a complex task.

    This function analyzes the complex task and generates a well-structured plan
    that outlines the necessary steps in the correct order to achieve the desired outcome.
    """
    plan_creation_prompt = (
        f"You are an expert in project planning, especially in the domains of bioinformatics, machine learning, and software development. "
        f"The user has provided a complex query that requires a step-by-step plan to execute. "
        f"Here is the query:\n\n"
        f"Query: '{query}'\n\n"
        "Please create a comprehensive sequential plan that outlines each necessary step clearly and logically. "
        "Ensure that the plan is detailed enough for implementation, and include any dependencies or prerequisites required for each step."
    )

    plan_response = chat_completion_api(plan_creation_prompt,system_prompt,tools=None)
    print("Hellllllll")
    return plan_response


def regenerate_plan(query, previous_plan):
    """
    Regenerates a detailed, step-by-step sequential plan based on user feedback and the previous plan.

    This function analyzes the user's new query in the context of the previously generated plan,
    and produces an updated, well-structured plan that incorporates the requested changes.
    """
    # Constructing the prompt to instruct the LLM to regenerate the plan
    plan_regeneration_prompt = (
        f"You are an expert in project planning, especially in the domains of bioinformatics, machine learning, and software development. "
        f"You have previously created a plan, but the user has requested changes. "
        f"Here is the original query and plan, followed by the user's new request:\n\n"
        f"Previous Plan:\n{previous_plan}\n\n"
        f"User's New Request: '{query}'\n\n"
        "Please regenerate the plan, incorporating the user's feedback and new requirements. "
        "Create a comprehensive sequential plan that outlines each necessary step clearly and logically. "
        "Ensure that the updated plan is detailed enough for implementation, and include any new or modified dependencies or prerequisites. "
        "Highlight the changes made from the previous plan."
    )

    # Making the LLM call to get the regenerated sequential plan
    updated_plan_response = chat_completion_api(plan_regeneration_prompt,system_prompt,tools=None)

    # Returning the detailed regenerated sequential plan
    return updated_plan_response






def calling():
    query1 = "Hi how are you?"
    query2 = "Can you write a simple for loop in python ?"
    query3 = "generate a plan for building a neural network ? also in the end append what function call you used from the set of tools ?"
    query4 = "I have uploaded 16 files, 4 files per cell-type. For each cell type there are three negative sequence files and one positive sequence file. Build a convolution neural network based model to classify positive and negative DNA sequences. For evaluation results, plot the area under precision recall curve and area under the receiver operator characteristic curve"
    system_prompt = """
    You are MasterAgent, an AI assistant specialized in bioinformatics. Your responsibilities include:

    1. User Interaction:
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


    4. Sequential Planning: -> use the function create_sequential_plan for this kind of queries please 
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


    5. Iterative Optimization:
        - Refine and optimize plans based on user feedback.

    6. Coder Agent Interaction:
        - Coordinate with the CoderAgent for tasks requiring extensive code generation or technical implementation.


    points : 

    1. Make plans only if its a big tasks or complex tasks not for smaller interactions and basic code.
    2. Use the supplied tools to assist the user.I will provide you with some tools which you need to use to assist the user
    3. MOST IMPORTANT : ANY LARGER QUERIES OR LARGER CODE GENERATION FIRST YOU NEED TO GENERATE A PLAN USING THE CREATE SEQUENTIAL PLAN FUNCTIONS
  
    """
    

    prev_plan = """
Step 1: Data Preparation\n- **1.1 Load Files**\n  - Use a programming language (e.g., Python) with libraries like `pandas` to load the 16 files containing DNA sequences.\n  \n- **1.2 Separate by Cell Type**\n  - Identify the cell types from the filenames or internal metadata and separate the sequences accordingly.\n  \n- **1.3 Label Sequences**\n  - Assign labels to the sequences: 'positive' for the sequences known to be associated with a condition and 'negative' for those that are not. \n  - This may involve creating a lookup table or a dictionary mapping filenames to labels if not embedded in the sequences.\n\n#### Step 2: Data Preprocessing\n- **2.1 Encode DNA Sequences**\n  - Implement one-hot encoding for the DNA sequences. For example, represent 'A' as [1, 0, 0, 0], 'C' as [0, 1, 0, 0], 'G' as [0, 0, 1, 0], and 'T' as [0, 0, 0, 1].\n  \n- **2.2 Split Data**\n  - Split the data into training, validation, and test sets. A typical split might be 70% training, 15% validation, and 15% testing.\n  - Ensure that the splits are stratified based on labels to maintain the ratio of positive and negative samples.\n\n#### Step 3: Model Design\n- **3.1 Choose Framework**\n  - Select a deep learning framework like TensorFlow or PyTorch for building the CNN model.\n\n- **3.2 Construct CNN Architecture**\n  - Design the CNN architecture, including:\n    - Input layer suitable for the one-hot encoded sequences.\n    - Convolutional layers (specify kernel sizes, number of filters).\n    - Activation functions (e.g., ReLU).\n    - Pooling layers (e.g., max pooling).\n    - Fully connected layers leading to the output layer.\n    - Output layer with activation function (e.g., sigmoid for binary classification).\n\n#### Step 4: Training\n- **4.1 Compile Model**\n  - Compile the model with appropriate loss functions (e.g., binary cross-entropy), optimizer (e.g., Adam), and metrics (e.g., accuracy).\n\n- **4.2 Training Process**\n  - Train the model using the training set and validate using the validation set. \n  - Set appropriate epochs and batch size and implement early stopping if necessary based on validation loss.\n\n- **4.3 Monitor Training**\n  - Capture training and validation loss/accuracy over epochs for later analysis.\n\n#### Step 5: Evaluation\n- **5.1 Evaluate on Test Set**\n  - Evaluate the trained model using the unseen test set to obtain overall performance metrics such as accuracy.\n\n- **5.2 Calculate Additional Metrics**\n  - Calculate precision and recall using the confusion matrix.\n  \n- **5.3 Generate ROC and Precision-Recall Curves**\n  - Use appropriate functions to plot the ROC curve and the Precision-Recall curve.\n\n#### Step 6: Visualization\n- **6.1 Plot and Interpret Curves**\n  - Plot the area under the Precision-Recall curve and ROC curve, using libraries such as `matplotlib` or `seaborn`.\n\n- **6.2 Summarize Results**\n  - Provide a summary of the model performance, discussing the implications of precision and recall, specifically in the context of DNA sequence classification.\n\n### Dependencies and Prerequisites\n- **Software Requirements:**\n  - Python environment with libraries: `numpy`, `pandas`, `tensorflow` or `pytorch`, `matplotlib`, `seaborn`, `scikit-learn`.\n  \n- **Data Availability:**\n  - Ensure access to the 16 DNA sequence files and all necessary metadata regarding cell types and labeling.

"""


    query5 = f"""
    query : So in the plan I want the training data to be 80% split and then validation should be 10% and then the test should be 10% . 
    change the plan accordingly please . 
    here is the previous plan : {prev_plan}
    """
    response =  process_query(query4,tools,system_prompt)
    print(response)

    if response.tool_calls:
    # Step 2: Extract the function name and arguments
        tool_call = response.tool_calls[0]  # Assuming there's only one tool call
        print("tool_call : ",tool_call)
        function_name = tool_call.function.name
        function_arguments = tool_call.function.arguments

        # Step 3: Convert the arguments from JSON to a Python dictionary
        args_dict = json.loads(function_arguments)
        print("args_dict : ",args_dict)

        # Step 4: Call the appropriate function
        # Assuming the functions are defined in the current namespace
        if function_name in globals():
            result = globals()[function_name](**args_dict)
            print("result :  " , result)
        else:
            print(f"Function {function_name} not found.")
    else:
        print("No tool call in the response.")


calling()