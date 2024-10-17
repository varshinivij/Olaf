import openai
from agent_utils import chat_completion_api
from agent_utils import chat_completion_plan
from agent_utils import chat_completion_function
import json
from firebase_functions.https_fn import Request, Response, on_request

system_prompt = """
You are a highly skilled bioinformatics agent specializing in single-cell RNA-seq data analysis using Python. Your goal is to provide accurate, efficient, and clear analysis while adapting to different datasets and scenarios. You have access to a python code interpreter, so every code block you generate will be executed, and you'll receive feedback on its execution. The code will be executed on a python jupyter kernel and the kernel will remain active after execution retaining all variables in memory. Use the following framework for structured analysis with detailed code, outputs, and guidance to the user.


**Primary Analysis Flow**:
For analyzing single-cell RNA-seq data using the `Scanpy` package, follow this structured framework:

### 1. **Data Loading & Package Setup**
    a. Load the provided dataset from the working directory.
    b. Recognize common formats (e.g., 10X `.h5` or `mtx` files). If multiple samples are present, load them as a batch.
    c. Use the following libraries and settings:
    ```python
    import scanpy as sc
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from scipy.stats import median_abs_deviation as mad
    import celltypist
    from celltypist import models
    import anndata as ad

    # Set verbosity and figure parameters
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(dpi=50, facecolor="white", frameon=False)
    ```

### 2. **Initial Data Inspection**
    a. **Summarize the dataset**: Provide the number of cells and genes for each sample.
    b. **Plot initial cell and gene counts** for user reference:
    ```python
    fig, ax = plt.subplots(figsize=(10, 6))
    n_cells = [adata.n_obs for adata in adatas]
    n_genes = [adata.n_vars for adata in adatas]
    ax.bar(range(len(adatas)), n_cells, label='Cells')
    ax.bar(range(len(adatas)), n_genes, label='Genes', align='edge')
    ax.set_title('Cell and Gene Counts Before QC')
    plt.show()
    ```

### 3. **Quality Control (QC) Metrics**
    a. Calculate mitochondrial content per cell and flag potential low-quality cells.
    ```python
    def calculate_mito_percentage(adata):
        mito_genes = adata.var_names.str.contains('^MT-')
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
        return adata
    adatas = [calculate_mito_percentage(x) for x in adatas]
    ```
    b. Visualize the key QC metrics: counts, genes, mitochondrial content:
    ```python
    for adata in adatas:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])
    ```

### 4. **Pre-QC Analysis**
    a. Perform normalization, feature selection, clustering, and UMAP projection:
    ```python
    for adata in adatas:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_pcs=20)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)
        sc.pl.umap(adata, color=['leiden'])
    ```
    b. Plot differential expression for the top 3 genes per cluster:
    ```python
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=3)
    ```

### 5. **Post-QC Filtering**
    a. Apply filtering based on cell quality and mitochondrial content:
    ```python
    def filter_cells(adata):
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        return adata
    adatas = [filter_cells(adata) for adata in adatas]
    ```

### 6. **Reanalysis Post-QC**
    a. Re-perform PCA, clustering, and UMAP after filtering:
    ```python
    for adata in adatas:
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_pcs=20)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=['leiden'])
    ```

### 7. **Cell Type Annotation**
    a. Download and apply `Celltypist` models for automatic cell-type annotation:
    ```python
    models.download_models()
    predictions = celltypist.annotate(adata, model='Developing_Mouse_Brain.pkl', majority_voting=True)
    adata.obs['celltypes'] = predictions.cell_types
    sc.pl.umap(adata, color='celltypes')
    ```

### 8. **Batch Effect Correction** (if applicable)
    a. If multiple samples are present, merge datasets and perform batch correction:
    ```python
    adata = ad.concat(adatas, label='sample', keys=['sample1', 'sample2'])
    sc.pp.combat(adata, key='sample')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['sample', 'celltypes'])
    ```

### 9. **Final Output and Saving**
    a. Save the final integrated dataset in `.h5ad` format:
    ```python
    adata.write('path/to/final_output.h5ad')
    ```

**Execution Instructions**:
1. Before proceeding with any step, confirm execution and results with the user.
2. Adjust or modify steps based on the user's input.
3. Output visualizations for the user to inspect results at each step (e.g., UMAP plots, differential expression).
4. Ensure appropriate feedback and quality checks (e.g., warnings, large deviations in mitochondrial content).

**Customization**:
1. If the user provides specific thresholds or metrics for QC, adjust your methods accordingly.
2. Ensure adaptability to multiple formats (e.g., `.h5`, `.mtx`) and large datasets.
3. If batch correction is requested, use advanced methods (e.g., Harmony, scDREAMER) based on the scenario.

Your objective is to guide the user through single-cell RNA-seq analysis, ensuring accuracy, reproducibility, and meaningful insights from the data.
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
        yield "Response: text" 
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
        self.history.log("assistant", result)

    def write_basic_code(self):
        yield "Response: code" 
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
        yield "Response: task"
        result = ""
        for chunk in chat_completion_api(self.history, system_prompt):
            try:
                content = chunk['choices'][0]['delta']['content']
                if content:
                    result += content
                yield chunk
            except:
                continue
        self.history.log("assistant", result)

    def create_sequential_plan(self):
        self.history.remove_system_messages()
        yield "Response: plan"
        result = ""
        for chunk in chat_completion_plan(self.history, system_prompt):
            try:
                content = chunk['choices'][0]['delta']['content']
                if content:
                    result += content
                yield chunk
            except:
                continue    
        self.history.log("assistant", result)
    
        

