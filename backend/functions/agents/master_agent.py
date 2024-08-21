from history import History
from agent_utils import chat_completion

system_prompt = """
You are `MasterAgent`, an experienced planner and requirements generator with expertise in bioinformatics and neural network modeling. Your task is to generate detailed, structured plans based on the user's query. Focus on the following aspects:
- Class structures, including relationships and components.
- Neural network architecture, layers, and connections.
- Preprocessing steps specific to genomic data.
- Essential methods/functions for data loading, preprocessing, and management.
- Integration of all components for a cohesive workflow.
"""

system = {
        "role": "system",
        "content": system_prompt
      }

class MasterAgent:

    def __init__(self):
        pass
        
    def plan(self, history):
        """
        This method takes a query, sends it to the LLM instance, and retrieves a structured plan or requirements.
        """
        hidden_prompt = self.generate_prompt(history.get_history())
        history.log("user", hidden_prompt)
        response = chat_completion(history)
        return response
    
    def re_plan(self, query, previous_plan):
        """
        This method takes a query and the previous plan, sends them to the LLM instance,
        and retrieves a new structured plan or requirements.
        """
        prompt = f"The user was not satisfied with the following plan:\n{previous_plan}\n\nPlease generate an improved plan based on the original query:\n{query}"
        response = chat_completion(self.history)
        return response

    def generate_prompt(self, query):
        """
        Generates a one-shot prompt for the LLM based on the user's query.
        """
        one_shot_example = (
            "Example Query: 'Generate detailed requirements for a neural network with LSTM for genomic data processing.'\n"
            "Example Response:\n"
            "{\n"
            "  'Class': 'GenomicDataProcessor',\n"
            "  'Class Name': 'GenomicDataProcessor',\n"
            "  'Language': 'Python',\n"
            "  'Additional Information': 'The class should include methods to load, preprocess, and split genomic data suitable for input into a neural network.',\n"
            "  'Components': [\n"
            "    {'Name': 'DataLoader', 'Parameters': ['data_path', 'batch_size']},\n"
            "    {'Name': 'Preprocessor', 'Parameters': ['normalization', 'data_augmentation']}\n"
            "  ],\n"
            "  'Methods': [\n"
            "    {'Name': 'load_data', 'Description': 'Loads data from the specified path.', 'Inputs': ['data_path'], 'Outputs': ['data']},\n"
            "    {'Name': 'preprocess_data', 'Description': 'Preprocesses the data for training.', 'Inputs': ['data'], 'Outputs': ['processed_data']}\n"
            "  ]\n"
            "}\n\n"
            f"Generate a detailed plan and requirements for the following query:\n\n{query}\n\n"
            "Include comprehensive information on the class structure, components, methods, and any other relevant details. Ensure that the requirements are precise, actionable, and aligned with best practices for bioinformatics and neural network modeling."
        )
        return one_shot_example