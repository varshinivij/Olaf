import requests
import re
import json

class MasterAgent:
    def __init__(self, model="GPT"):
        self.task_agents = 0
        self.task_count = 0
        self.model = model
        self.system_message = """You are a master bioinformatician agent responsible for planning and executing analysis on single-cell data. You have access to a code interpreter where every code you generate will be executed, and you will receive the output. Your goal is to ensure accurate and efficient analysis. Utilize the following framework:

        1. First, create a step-by-step plan as a numbered list to break down the problem into individual tasks that you can perform. The following is a skeleton for analyzing single-cell data:

        (a) Load the file in the working directory based on the file extension and provide a summary of the contents within the file.
        (b) If it is an h5ad file, provide the number of cells within each sample.
        (c) Generate and plot QC metrics.
        (d) Perform QC filtering with the default cutoffs.
        (e) Normalize, scale, and plot the UMAP.

        2. Once you create the plan, execute each task one at a time. Wait for the code for each step to be executed and the user's confirmation before proceeding to the next step.

        Guidelines for specific steps:

       1. Before generating QC metrics, check gene names to identify MT genes. Use the following code as a framework:
       adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        sc.settings.set_figure_params(dpi=150, dpi_save=300, figsize=(12,3), facecolor='white')
        sc.pl.violin(adata, ['total_counts'], inner='box', size=0, groupby='SampleID', multi_panel=False, rotation=45)

      Similarly, plot n_genes_by_counts and pct_counts_mt when asked for QC.

      2. Filtering scRNA-seq data should involve retaining cells with at least 200 genes, excluding genes detected in fewer than 3 cells, and removing cells with total counts below 250 or above 30,000. Additionally, filter out cells with more than 6,000 genes and retain those with less than or equal to 2% mitochondrial gene counts to ensure high-quality data.

      3. Plotting UMAP should involve normalizing based on total counts, log-normalizing, identifying highly variable genes, scaling, computing PCA, neighborhood graph, UMAP, performing Leiden clustering, and plotting the UMAP.
      
      """

    def chat_completion(self, history):
        api_key = "REMOVED"
        url = "https://api.openai.com/v1/chat/completions"

        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }

        payload = {
            "model": "gpt-4o",
            "messages": history,
            "temperature": 0.1
        } 

        response = requests.post(url, headers=headers, json=payload)

        if response.status_code == 200:
            result = response.json()
            response_message = result['choices'][0]['message']['content']
            return response_message

        else:
            print(f"Error: {response.status_code}, {response.text}")
            return None
        
    def parse_output(self, output):
        code_block_pattern = r"```(.*?)```"
        code_matches = re.finditer(code_block_pattern, output, re.DOTALL)
        data = []
        last_end = 0
        for match in code_matches:
            start, end = match.span()
            text_content = output[last_end:start].strip()
            if text_content:
                data.append({"type": "text", "role": "assistant", "content": text_content})
            
            code_content = match.group(1).strip()
            if code_content.startswith('python'):
                code_content = code_content[6:].strip()
            data.append({"type": "code", "role": "assistant", "content": code_content})
            last_end = end

        remaining_text = output[last_end:].strip()
        if remaining_text:
            data.append({"type": "text", "role": "assistant", "content": remaining_text})

        return json.dumps(data, indent=4)