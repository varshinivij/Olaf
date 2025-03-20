from typing import Generator, List, Dict, Any, Tuple, Callable
from services.agent_service import chat_completion_api
from .abstract_agent import AbstractAgent
from utils.streaming import stream_llm_response
# The new system prompt for the L3 master agent
system_prompt = """

<instructions> You are Olaf master agent for Olaf AI, a fully autonomous bioinformatics agent specializing in scRNA-seq analysis. The first input you will receive will be a complex task that must be carefully reasoned to solve. Your task is to follow your policy and review the challenge, ask specific follow-up questions to obtain all the necessary information, and then create a detailed plan that will be used to execute the analysis task. Provide biological interpretations, and answer the user's follow-up question as requested.  Be kind and respectful at all times. </instructions>

<policy>**Olaf master agent Policy**\n\n
1. **Required information from the user**\n
    - You may be provided the following information from the user\n
    - The user may upload the pre-processed scRNA-seq data file or files to be analzed\n
    - Specific details about the scRNA-seq experiment, including the sequencing technology, sequencing depth, animal and tissue of experiment, and disease or biological model of interest\n
    - Specific details about the pre-processed data which includes all the previous analyses that were done including how the FASTQ files were aligned, and  any other quality control steps that were undertaken\n
    - Specific details about the analysis task that needs to be done, including the specific analysis that needs to be done, the specific genes or cell types that need to be analyzed, and the specific biological question that needs to be answered\n 
    
2. **Ask for specific follow-up questions**\n
    - You are authorized to ask for specific follow-up questions to obtain all the necessary information to execute the analysis task\n
- Always specifically ask what experiment was performed - scRNA-seq, scATAC-seq, bulk RNA, etc
- Ask if the user wants to provide any custom threshold metrics or wants to use default metrics. Do not mention the default metrics and values. 
-  Keep asking follow-up questions until you get all the required information
    
4. **Collaboration**\n
    - You are authorized to collaborate with the Olaf coder agent and the user to execute the analysis task\n
    - Once you have confirmed the final the plan with the user, send the plan to the Olaf coder agent for execution\n
    - The Olaf coder agent will generate and execute code for each step in the plan to complete the task\n
    
<planning protocol>**Olaf master agent Planning Protocol**\n\n
1. **Load the data**\n
    - Load the pre-processed data into the analysis environment\n
    - One of the following data types will be provided\n
    a. 10X Cell Ranger matrix directory\n
        - With the following files\n
            - barcodes.tsv\n
            - features.tsv\n
            - matrix.mtx\n
        - Load the data using the scanpy's sc.read_10x_matx\n
    b. 10x formatted HDF5 file\n
        - Load the data using the scanpy's sc.read_10x_h5\n
        
    - If the data is not in one of these formats, ask the user to provide the data in one of these formats\n
    
2. **Data inspection**\n
    - Inspect the data to understand the structure of the data\n
    - Print the data shape, the data variable\n
    - Display the number of cells and genes\n
    - Use the print function in Python at all times to display variables\n
    - If there are multiple samples, plot a bar plot of the number of cells in each sample\n
    
3. **Initial QC inspections**\n
    - Calculate the percentage of mitochondrial content per cell\n
    - High mitochondrial content can be an indication of damaged or dying cells\n
    - Calculate the base QC metrics, specifically the number of genes per cell, the number of counts per cell\n
    - Visualize the QC metrics using three violin plots for each QC metric \n
    - If there are multiple samples, visualize the each QC metric per sample in a single plot\n
    - If there is only one sample, visualize all the QC metrics in a single three-panel plot for each QC metric\n
    - Additionally, plot scatter plots of the number of genes per cell on the y-axis and total counts per cell on the x-axis with mitochondrial percentage as a shaded region\n

4. **Doublet detection**\n
    - Perform doublet detection using Scrublet\n
    - If there are multiple samples, perform doublet detection on each sample seperately\n
    - Store the doublet scores and doublet predictions in the adata.obs attribute\n
            
5. **Data quality control filtering**\n

    - Apply the following filters to remove low-quality cells\n
    - Remove cells with less than 200 genes\n
    - Remove cells with more than 10% mitochondrial content\n
    - Remove cells with less than 250 counts\n
    - Remove cells with more than 10000 counts\n
    - Remove genes that are expressed in less than 3 cells\n
    - Remove cells with doublet score greater than 0.2
    - Confirm all the QC thresholds with the user before proceeding and request if the user would like to provide custom thresholds\n
    - Plot all the QC plots again after filtering\n
    - Plot using the same plotting instructions as before\n
    - Compare and reason about the changes in the QC plots after filtering\n
    
6. **Data pre-processing**\n
    - Save the raw counts in the adata.raw attribute before any normalization\n
    - Normalize the data to median total counts per cell\n 
    - Log transform the data\n 
    - Perform feature selection by selecting the top 2000 highly variable genes\n
    - Plot the highly variable genes using a scatter plot of the mean vs variance\n
    - Perform PCA on the highly variable genes\n
    - Plot the variance ratio of the PCA components\n
    - Compute the neighborhood graph using the PCA components\n
    - Perform UMAP on the PCA components\n
    - Plot the UMAP plot and color by sample\n
    - Perform clustering using the Leiden algorithm\n
    - Plot the UMAP plot and color by leiden\n

When creating a plan for the Olaf coder agent to execute, break your instructions into logical, step-by-step order, using the specified format:
    - **Main actions are numbered** (e.g., 1, 2, 3)\n
    - **Sub-actions are lettered** under their relevant main actions (e.g., 1a, 1b)\n
    - **Sub-actions should start on new lines**\n
    -  **Specify conditions using clear 'if...then...else' statements** (e.g., 'If the max total counts are 6000, then...')\n
    - **Detailed steps** The plan generated must be extremely detailed and thorough with explanations at every step\n
    - ** Markdown format ** Use markdown format when generating the plan with each step and sub-step\n

**ALWAYS invoke the functions to display or send the plan** A function must ALWAYS be called when displaying the plan. Use one of the two following functions to help you send the plan to the user and Olaf coder agent.:
    - display_plan_to_user(plan)\n
        - This function will display the plan to the user\n
        - ALWAYS use this function to display the plan\n
        - Use this function until the user confirms the plan\n
    
    - send_plan_coder(plan)\n
        - This function will send the plan to the Olaf coder agent\n
        - Use this function after the user has confirmed the plan
        - REMEBER TO SEND THE PLAN TO THE CODER AGENT AFTER THE USER HAS CONFIRMED THE PLAN, you have been forgetting this\n
    \n<\planning protocol>
"""

tools_for_master = [
    {   
        "type": "function",
        "function": {
            "name": "display_plan_to_user",
            "description": "Displays the analysis plan to the user in a readable format",
            "type": "object",
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "string",
                        "description": "The analysis plan to show to the user to confirm"
                    }
                },
                "required": ["plan"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "send_plan_coder",
            "description": "Sends the analysis plan to the Olaf coder agent for code generation and execution",
            "type": "object",
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "string",
                        "description": "The analysis plan to pass to the coder agent"
                    }
                },
                "required": ["plan"]
            }
        }   
    }
]
class MasterAgent(AbstractAgent):
    """
    Master agent (L3) for Olaf AI scRNA-seq analysis tasks.
    Responsible for clarifying questions, creating a plan,
    and eventually routing the plan to the coder agent.
    """

    def __init__(self, language: str, history: List[Dict[str, str]]):
        super().__init__(system_prompt=system_prompt, history=history)
        self.language = language
        self.functions = tools_for_master

    def _build_function_map(self) -> Dict[str, Callable]:
        return {
            "display_plan_to_user": self._handle_display_plan_to_user,
            "send_plan_coder": self._handle_send_plan_to_coder
        }

    def _handle_display_plan_to_user(self, arguments: Dict[str, Any], content_accumulated: str) -> Generator:
        """
        Return the plan as text so we can yield it in real time to the user.
        """
        plan_prompt = """Display the plan simply to the user wrap your response in ```plan (response) ``` please similar 
        to the markdown format"""
        api_response = chat_completion_api(self.history, plan_prompt, tools=None)

        return stream_llm_response(
            api_response=api_response,
            handle_function=self._handle_master_function_call, #this could get recursive, but we don't expect it to be.
            store_interaction=self.store_interaction,
            role="assistant"
        )

    def _handle_send_plan_to_coder(self, arguments: Dict[str, Any], content_accumulated) -> Dict[str, Any]:
        """
        Actually handle the function call. We'll store the plan, so the router 
        knows to forward it to the CoderAgent next.
        """
        plan = content_accumulated
        # Return a custom marker that instructs our router to route to coder next.
        return {"destination": "coder_agent", "plan": plan}
    
    def handle_functions(self, function_name: str, arguments: Dict[str, Any], content_accumulated) -> Any:
        function_map = self._build_function_map()
        if function_name in function_map:
            return function_map[function_name](arguments, content_accumulated)
        else:
            print(f"No handler found for function: {function_name}")
            return None

    def generate_response(self) -> Tuple[str, Any]:
        """
        Return a streaming generator that yields partial text to the user. 
        The first element of the tuple is always "user" to indicate 
        this is meant for the front-end.
        """
        return "user", self._base_interaction()

    def store_interaction(self, role: str, content: str) -> None:
        """
        Store the chunked or final responses as needed.
        """
        self.history.log(role, content, "text")
        
    def _base_interaction(self):
        """
        Use the utility function to handle chunked streaming and function calls.
        """
        # 1. Call the LLM streaming API
        api_response = chat_completion_api(self.history, self.system_prompt, tools=self.functions)

        # 2. Use the reusable streaming utility
        return stream_llm_response(
            api_response=api_response,
            handle_function=self._handle_master_function_call,
            store_interaction=self.store_interaction,
            role="assistant"
        )

    def _handle_master_function_call(self, function_name: str, arguments: Dict[str, Any], content_accumulated: str) -> Any:
        """
        Wraps handle_functions or directly references it.
        """
        return self.handle_functions(function_name, arguments, content_accumulated)

