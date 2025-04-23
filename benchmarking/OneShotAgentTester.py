import argparse
import os
import sys
import json
import re
import shlex
import time
from pathlib import Path
import subprocess # Added for docker cp

# --- Dependency Imports ---
try:
    from dotenv import load_dotenv
except ImportError:
    print("Error: python-dotenv library not found. Please install it: pip install python-dotenv", file=sys.stderr)
    sys.exit(1)

try:
    from openai import OpenAI, APIError
except ImportError:
    print("Error: openai library not found. Please install it: pip install openai", file=sys.stderr)
    sys.exit(1)

try:
    # Assumes benchmarking_sandbox_management.py is in a 'sandbox' subdirectory
    sandbox_dir = os.path.join(os.path.dirname(__file__), 'sandbox')
    sys.path.insert(0, sandbox_dir)
    # Import both the class and the constant
    from benchmarking_sandbox_management import SandboxManager, CONTAINER_NAME as SANDBOX_CONTAINER_NAME
except ImportError as e:
    print(f"Error: Could not import SandboxManager or CONTAINER_NAME from {sandbox_dir}.", file=sys.stderr)
    print("Ensure benchmarking_sandbox_management.py is present in the 'sandbox' directory.", file=sys.stderr)
    print(f"Details: {e}", file=sys.stderr)
    sys.exit(1)
finally:
    # Clean up sys.path modification if SandboxManager was imported
    if 'benchmarking_sandbox_management' in sys.modules:
        sys.path.pop(0)


# Optional: Use rich for better formatting
try:
    from rich.console import Console
    from rich.prompt import Prompt, Confirm
    from rich.panel import Panel
    from rich.syntax import Syntax
    from rich.table import Table
    HAS_RICH = True
except ImportError:
    HAS_RICH = False
    # Simple print/input fallback if rich is not installed
    class Console:
        def print(self, *args, **kwargs): print(*args)
    class Prompt:
        @staticmethod
        def ask(prompt, choices=None, default=None):
            p_text = f"{prompt} "
            if choices: choices_str = '/'.join(choices); p_text += f"({choices_str}) "
            if default: p_text += f"[{default}] "
            return input(p_text).strip()
        @staticmethod
        def get_input(console, prompt, password=False): # Match rich signature somewhat
             return input(f"{prompt}: ")
    class Confirm:
         @staticmethod
         def ask(prompt, default=False):
             val = input(f"{prompt} [y/N] " if not default else f"{prompt} [Y/n] ").lower().strip()
             if not val: return default
             return val == 'y'
    class Panel: # Dummy Panel
        def __init__(self, content, title="", border_style=""): self.content=content; self.title=title
        def __rich_console__(self, console, options): yield self.title; yield self.content
    class Syntax: # Dummy Syntax
        def __init__(self, code, lexer, theme="", line_numbers=False): self.code = code; self.lexer = lexer
        def __rich_console__(self, console, options): yield f"--- Code ({self.lexer}) ---\n{self.code}\n--- End Code ---"
    class Table: # Basic fallback Table class
        def __init__(self, title=""): self._title=title; self._rows=[]; self._columns=[]
        def add_column(self, header, style="", justify="left", no_wrap=False): self._columns.append(header)
        def add_row(self, *items):
             if len(items) != len(self._columns): raise ValueError("Row items != columns")
             self._rows.append(items)
        def __rich_console__(self, console, options):
             yield self._title;
             if self._columns:
                 yield "\t".join(self._columns)
                 for row in self._rows: yield "\t".join(map(str, row))
        def print_table(self, console): # Custom print method if rich not available
             console.print(self._title)
             if self._columns:
                 col_widths = [len(h) for h in self._columns]
                 for row in self._rows:
                     for i, item in enumerate(row): col_widths[i] = max(col_widths[i], len(str(item)))
                 header = "  ".join(f"{h:<{w}}" for h, w in zip(self._columns, col_widths))
                 separator = "-" * len(header)
                 console.print(header); console.print(separator)
                 for row in self._rows:
                     console.print("  ".join(f"{str(item):<{w}}" for item, w in zip(row, col_widths)))

# --- Constants ---
SCRIPT_DIR = Path(__file__).parent.resolve()
DATASETS_DIR = SCRIPT_DIR / "datasets"
ENV_FILE = SCRIPT_DIR / ".env"
SANDBOX_DATA_PATH = "/home/user/data.h5ad" # Where data will be copied inside container

# --- Configuration Loading ---
console = Console()
load_dotenv(dotenv_path=ENV_FILE)
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

if not OPENAI_API_KEY:
    console.print(f"[bold red]Error:[/bold red] OPENAI_API_KEY not found in {ENV_FILE}.")
    console.print("Please run the 'make_benchmarking_env.sh' script first.")
    sys.exit(1)

try:
    openai_client = OpenAI(api_key=OPENAI_API_KEY)
    # Test connection (optional, but good practice)
    # openai_client.models.list()
    # console.print("OpenAI client initialized successfully.")
except Exception as e:
    console.print(f"[bold red]Error initializing OpenAI client:[/bold red] {e}")
    sys.exit(1)


# --- Helper Functions ---
def extract_python_code(text):
    """Extracts the first Python code block from text."""
    match = re.search(r"```python\s*([\s\S]+?)\s*```", text, re.MULTILINE)
    if match:
        return match.group(1).strip()
    return None

def display_message(role, content):
    """Displays messages with nice formatting."""
    # Existing display logic... (unchanged)
    if role == "system":
        console.print(Panel(content, title="SYSTEM PROMPT", border_style="dim blue"))
    elif role == "user":
        # Check if it's the special code execution result message
        if content.startswith("Code execution result:\n"):
             console.print(Panel(content, title="CODE EXECUTION RESULT (Sent as User)", border_style="yellow"))
        else:
             console.print(Panel(content, title="USER (Input Prompt)", border_style="blue"))
    elif role == "assistant":
        code = extract_python_code(content)
        if code:
            text_part = re.sub(r"```python\s*([\s\S]+?)\s*```", "", content, count=1).strip()
            if text_part:
                 console.print(Panel(text_part, title="ASSISTANT (Text)", border_style="green"))
            if HAS_RICH:
                 console.print(Panel(Syntax(code, "python", theme="default", line_numbers=True), title="ASSISTANT (Code)", border_style="green"))
            else:
                 console.print(f"--- ASSISTANT (Code) ---\n{code}\n--- End Code ---")
        else:
            console.print(Panel(content, title="ASSISTANT (Text Only)", border_style="green"))
    # Remove the explicit 'tool' role display as we are sending it as 'user'
    # elif role == "tool":
    #      console.print(Panel(content, title="CODE EXECUTION RESULT", border_style="yellow"))
    else:
        console.print(f"[bold]{role.upper()}:[/bold]\n{content}")
    console.print("-" * 20) # Separator

# --- Core Logic Functions ---

def get_agent_prompts():
    """Gets agent prompt(s) based on user input method."""
    prompts = {}
    while True:
        console.print("\n[bold cyan]Select Agent Prompt Input Method:[/bold cyan]")
        console.print("  1. Paste prompt directly into the terminal.")
        console.print("  2. Provide path to a single .txt file.")
        console.print("  3. Provide path to a folder containing .txt prompt files.")
        choice = Prompt.ask("Enter choice (1/2/3)", choices=["1", "2", "3"], default="1")
        if choice == "1":
            console.print("Paste your prompt below. Press Ctrl+D (Unix) or Ctrl+Z+Enter (Windows) when done:")
            try:
                prompt_text = sys.stdin.read().strip()
                if prompt_text: prompts["pasted_prompt"] = prompt_text; return prompts
                else: console.print("[yellow]No prompt pasted. Please try again.[/yellow]")
            except EOFError: console.print("\n[yellow]No prompt pasted. Please try again.[/yellow]")
        elif choice == "2":
            file_path_str = Prompt.ask("Enter the path to the .txt prompt file")
            file_path = Path(file_path_str).resolve()
            if file_path.is_file() and file_path.suffix.lower() == ".txt":
                try:
                    prompt_text = file_path.read_text(encoding='utf-8').strip()
                    if prompt_text: prompts[file_path.stem] = prompt_text; return prompts
                    else: console.print(f"[yellow]File '{file_path}' is empty.[/yellow]")
                except Exception as e: console.print(f"[red]Error reading file '{file_path}': {e}[/red]")
            else: console.print(f"[red]Invalid path or not a .txt file: '{file_path_str}'[/red]")
        elif choice == "3":
            folder_path_str = Prompt.ask("Enter the path to the folder containing .txt prompt files")
            folder_path = Path(folder_path_str).resolve()
            if folder_path.is_dir():
                txt_files = list(folder_path.glob("*.txt"))
                if not txt_files: console.print(f"[yellow]No .txt files found in folder '{folder_path_str}'.[/yellow]"); continue
                for file_path in txt_files:
                    try:
                        prompt_text = file_path.read_text(encoding='utf-8').strip()
                        if prompt_text: prompts[file_path.stem] = prompt_text
                        else: console.print(f"[yellow]Skipping empty file: '{file_path.name}'[/yellow]")
                    except Exception as e: console.print(f"[red]Error reading file '{file_path.name}': {e}[/red]")
                if prompts: console.print(f"Found {len(prompts)} non-empty prompt files."); return prompts
                else: console.print("[yellow]No valid, non-empty prompts found in the folder.[/yellow]")
            else: console.print(f"[red]Invalid path or not a directory: '{folder_path_str}'[/red]")

def select_dataset():
    """Scans datasets directory and prompts user for selection."""
    if not DATASETS_DIR.is_dir():
        console.print(f"[bold red]Error:[/bold red] Datasets directory not found at '{DATASETS_DIR}'")
        console.print("Please ensure datasets are downloaded using 'czi_browser.py download ...'")
        return None, None
    datasets = []
    for h5ad_path in DATASETS_DIR.glob("*.h5ad"):
        json_path = h5ad_path.with_suffix(".json")
        if json_path.is_file():
            try:
                with open(json_path, 'r', encoding='utf-8') as f: metadata = json.load(f)
                datasets.append({ "h5ad_path": h5ad_path, "json_path": json_path, "metadata": metadata, "display_name": metadata.get("dataset_title", h5ad_path.stem)})
            except Exception as e: console.print(f"[yellow]Warning: Could not load metadata for '{h5ad_path.name}': {e}[/yellow]")
        else: console.print(f"[yellow]Warning: Missing metadata file for '{h5ad_path.name}'. Skipping.[/yellow]")
    if not datasets:
        console.print(f"[bold red]Error:[/bold red] No valid datasets found in '{DATASETS_DIR}'")
        return None, None
    console.print("\n[bold cyan]Available Datasets:[/bold cyan]")
    table = Table(title="Select a Dataset")
    table.add_column("Index", style="dim", justify="right")
    table.add_column("Dataset Title / Filename", style="green")
    table.add_column("Cell Count", style="magenta", justify="right")
    table.add_column("Organism", style="blue")
    for i, ds in enumerate(datasets):
        meta = ds["metadata"]; cell_count = meta.get('cell_count', 'N/A'); organism = ", ".join(meta.get('organism', [])) if isinstance(meta.get('organism'), list) else meta.get('organism', 'N/A')
        try: cell_count_str = f"{int(cell_count):,}" if cell_count != 'N/A' else 'N/A'
        except (ValueError, TypeError): cell_count_str = str(cell_count)
        table.add_row(str(i + 1), ds["display_name"], cell_count_str, organism)
    if HAS_RICH: console.print(table)
    else: table.print_table(console)
    while True:
        choice_str = Prompt.ask(f"Enter the index of the dataset to use (1-{len(datasets)})")
        try:
            choice_idx = int(choice_str) - 1
            if 0 <= choice_idx < len(datasets):
                selected_ds = datasets[choice_idx]
                console.print(f"Selected dataset: [green]{selected_ds['display_name']}[/green]")
                return selected_ds["h5ad_path"], selected_ds["metadata"]
            else: console.print(f"[red]Invalid index. Please enter a number between 1 and {len(datasets)}.[/red]")
        except ValueError: console.print("[red]Invalid input. Please enter a number.[/red]")

def get_code_tries():
    """Prompts user for the number of code execution attempts."""
    while True:
        tries_str = Prompt.ask("Enter the maximum number of code execution attempts for the agent", default="5")
        try:
            tries = int(tries_str)
            if tries > 0: return tries
            else: console.print("[red]Number of tries must be positive.[/red]")
        except ValueError: console.print("[red]Invalid input. Please enter an integer.[/red]")

def run_agent_test(agent_prompt_id, agent_prompt, dataset_h5ad_path, dataset_metadata, max_code_tries):
    """Runs a single agent test loop."""
    console.print(f"\n[bold cyan]----- Starting Test: '{agent_prompt_id}' ----- [/bold cyan]")
    console.print(f"Dataset: [green]{dataset_metadata.get('dataset_title', dataset_h5ad_path.stem)}[/green]")
    console.print(f"Max Code Tries: [yellow]{max_code_tries}[/yellow]")

    sandbox_manager = None
    conversation_history = []
    code_tries_left = max_code_tries

    try:
        # 1. Initialize and Start Sandbox
        console.print("\nInitializing Sandbox Manager...")
        sandbox_manager = SandboxManager()
        console.print("Starting sandbox container...")
        container = sandbox_manager.start_container()
        if not container:
            console.print("[bold red]Failed to start sandbox container. Aborting test.[/bold red]")
            return None

        # 2. Copy Dataset to Sandbox
        console.print(f"Copying dataset '{dataset_h5ad_path.name}' to sandbox ({SANDBOX_DATA_PATH})...")
        copy_command = ['docker', 'cp', str(dataset_h5ad_path), f"{SANDBOX_CONTAINER_NAME}:{SANDBOX_DATA_PATH}"]
        try:
            subprocess.run(copy_command, check=True, capture_output=True, text=True)
            console.print("[green]Dataset copied successfully.[/green]")
        except subprocess.CalledProcessError as e:
            console.print(f"[bold red]Error copying dataset to container:[/bold red]")
            console.print(f"Command: {' '.join(e.cmd)}")
            console.print(f"Stderr: {e.stderr}")
            raise

        # 3. Prepare Initial Agent Message
        system_message_content = f"""You are an AI assistant tasked with analyzing a single-cell transcriptomics dataset.
Your goal is to characterize this dataset based on its metadata and by executing Python code.
You have access to a Python environment within a sandbox. Standard libraries like pandas, numpy, scipy, scikit-learn, and anndata should be available.
The dataset is loaded into the sandbox at: {SANDBOX_DATA_PATH}
You can load it using anndata.read_h5ad('{SANDBOX_DATA_PATH}').

Dataset Metadata:
{json.dumps(dataset_metadata, indent=2)}

You have a maximum of {max_code_tries} attempts to execute Python code blocks.
When you want to execute code, enclose it in triple backticks with the language specified as python, like this:
```python
import anndata as ad
adata = ad.read_h5ad('{SANDBOX_DATA_PATH}')
# Your analysis code here
print(adata.shape)
```
I will run the code you provide and return the output (stdout and stderr). Use the output to inform your next step.
Focus on providing meaningful characterizations and insights based on the data and metadata. Plan your {max_code_tries} code executions wisely. Start by loading the data and examining its basic properties.
"""
        user_message_content = agent_prompt

        conversation_history = [
            {"role": "system", "content": system_message_content},
            {"role": "user", "content": user_message_content}
        ]
        display_message("system", system_message_content)
        display_message("user", user_message_content)

        # 4. Agent Interaction Loop
        while code_tries_left > 0:
            console.print(f"\n[bold]Sending request to OpenAI... (Code tries left: {code_tries_left})[/bold]")
            try:
                response = openai_client.chat.completions.create(
                    model="gpt-4o", # Or your preferred model
                    messages=conversation_history,
                    temperature=0.7,
                )
                assistant_message = response.choices[0].message
                assistant_content = assistant_message.content

                # Append assistant's response BEFORE processing code
                conversation_history.append({"role": "assistant", "content": assistant_content})
                display_message("assistant", assistant_content)

                # 5. Check for and Execute Code
                code_to_run = extract_python_code(assistant_content)
                if code_to_run:
                    console.print(f"\n[bold cyan]Executing Code (Attempt {max_code_tries - code_tries_left + 1}/{max_code_tries})...[/bold cyan]")
                    execution_output = sandbox_manager.run_code(code_to_run)
                    code_tries_left -= 1

                    # Prepare result message for history and display
                    # **MODIFICATION:** Send result back as 'user' role to avoid API error
                    user_feedback_content = f"Code execution result:\n"
                    if execution_output is not None:
                         # Limit output length sent back to OpenAI if necessary
                         max_output_len = 2000
                         if len(execution_output) > max_output_len:
                              user_feedback_content += f"--- STDOUT (Truncated) ---\n{execution_output[:max_output_len]}...\n--------------"
                         else:
                              user_feedback_content += f"--- STDOUT ---\n{execution_output}\n--------------"
                    else:
                         user_feedback_content += "[No standard output captured]"
                    # Note: stderr is printed by run_code but not added to history here.
                    # Could add stderr to user_feedback_content if needed.

                    # Append the execution result as a user message
                    conversation_history.append({"role": "user", "content": user_feedback_content})
                    # Display this feedback message (using the modified display_message)
                    display_message("user", user_feedback_content)

                    if code_tries_left == 0:
                        console.print("[bold yellow]Maximum code execution attempts reached.[/bold yellow]")
                        break # Exit loop

                else:
                    console.print("[yellow]No code block found in assistant's response this turn.[/yellow]")

            except APIError as e:
                console.print(f"[bold red]OpenAI API Error:[/bold red] {e}")
                # Attempt to print more details from the error object if available
                if hasattr(e, 'body') and e.body:
                     console.print(f"Error Body: {e.body}")
                break # Stop test on API error
            except Exception as e:
                console.print(f"[bold red]Error during agent interaction:[/bold red] {e}")
                break # Stop test on other errors

        console.print(f"\n[bold cyan]----- Test Finished: '{agent_prompt_id}' ----- [/bold cyan]")
        return conversation_history

    except Exception as e:
        console.print(f"[bold red]An error occurred during test setup or execution for '{agent_prompt_id}':[/bold red] {e}")
        return None
    finally:
        # 6. Stop and Cleanup Sandbox
        if sandbox_manager:
            console.print("\nStopping sandbox container...")
            if not sandbox_manager.stop_container():
                 console.print("[yellow]Warning: Could not cleanly stop/remove sandbox container.[/yellow]")



def main():
    parser = argparse.ArgumentParser(description="Run AI agent benchmarks against datasets in a sandbox.")
    args = parser.parse_args()
    console.print("[bold blue]Welcome to the One-Shot Agent Tester![/bold blue]")
    agent_prompts = get_agent_prompts()
    if not agent_prompts: console.print("[red]No agent prompts provided. Exiting.[/red]"); sys.exit(1)
    dataset_h5ad_path, dataset_metadata = select_dataset()
    if not dataset_h5ad_path or not dataset_metadata: console.print("[red]No dataset selected or available. Exiting.[/red]"); sys.exit(1)
    max_code_tries = get_code_tries()
    results = {}
    for prompt_id, prompt_text in agent_prompts.items():
        test_result = run_agent_test(prompt_id, prompt_text, dataset_h5ad_path, dataset_metadata, max_code_tries)
        results[prompt_id] = test_result
        if len(agent_prompts) > 1:
             if not Confirm.ask(f"\nTest for '{prompt_id}' finished. Continue with the next agent prompt?", default=True):
                  console.print("[yellow]Aborting remaining tests.[/yellow]"); break
             console.print("\n" + "="*40 + "\n"); time.sleep(1)
    console.print("\n[bold blue]All specified agent tests have concluded.[/bold blue]")
    # TODO: Process/save the 'results' dictionary

if __name__ == "__main__":
    main()
