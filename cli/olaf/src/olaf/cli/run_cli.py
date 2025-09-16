Got it — here’s your **final, conflict-free version** of `run_cli.py` with all 3 conflicts resolved and the correct imports:

```python
# olaf/cli/run_cli.py
import os
import re
import textwrap
from pathlib import Path
from typing import List, Tuple, cast, Optional
import subprocess
import json
from datetime import datetime

import typer
from rich.console import Console
from rich.prompt import Prompt, IntPrompt
from dotenv import load_dotenv
from olaf.config import DEFAULT_AGENT_DIR, ENV_FILE, OLAF_HOME

PACKAGE_ROOT = Path(__file__).resolve().parent.parent
PACKAGE_AGENTS_DIR = PACKAGE_ROOT / "agents"
PACKAGE_DATASETS_DIR = PACKAGE_ROOT / "datasets"
PACKAGE_AUTO_METRICS_DIR = PACKAGE_ROOT / "auto_metrics"

SANDBOX_DATA_PATH = "/workspace/dataset.h5ad"
SANDBOX_REF_DATA_PATH = "/workspace/reference.h5ad"

def _prompt_for_file(
    console: Console, user_dir: Path, package_dir: Path, extension: str, prompt_title: str
) -> Path:
    """
    Generic helper to find files in both user and package directories and prompt for a selection.
    """
    console.print(f"[bold]Select {prompt_title}:[/bold]")
    found_files = []
    seen_filenames = set()
    if user_dir.exists():
        for file_path in sorted(list(user_dir.glob(f"**/*{extension}"))):
            if file_path.name not in seen_filenames:
                found_files.append({"path": file_path, "label": "User"})
                seen_filenames.add(file_path.name)
    if package_dir.exists():
        for file_path in sorted(list(package_dir.glob(f"**/*{extension}"))):
            if file_path.name not in seen_filenames:
                found_files.append({"path": file_path, "label": "Package"})
                seen_filenames.add(file_path.name)
    if not found_files:
        console.print(f"[bold red]No '{extension}' files found.[/bold red]")
        raise typer.Exit(1)
    for i, file_info in enumerate(found_files, 1):
        console.print(f"  [cyan]{i}[/cyan]: {file_info['path'].name} [yellow]({file_info['label']})[/yellow]")
    choice_str = Prompt.ask("Enter the number of your choice", choices=[str(i) for i in range(1, len(found_files) + 1)])
    return found_files[int(choice_str) - 1]['path']

def _prompt_for_driver(console: Console, system: 'AgentSystem') -> str:
    """Prompts the user to select a driver agent from the loaded system."""
    console.print("[bold]Select a driver agent:[/bold]")
    agents = list(system.agents.keys())
    return Prompt.ask("Enter the name of the driver agent", choices=agents, default=agents[0])

def _prompt_for_benchmark_module(console: Console) -> Optional[Path]:
    """Finds and prompts the user to select an auto metric script."""
    console.print("[bold]Select a benchmark module (optional):[/bold]")
    
    modules = [
        m for m in PACKAGE_AUTO_METRICS_DIR.glob("*.py")
        if m.name not in ["__init__.py", "AutoMetric.py"]
    ]
    
    if not modules:
        console.print("[yellow]No benchmark modules found.[/yellow]")
        return None
        
    for i, mod in enumerate(modules, 1):
        console.print(f"  [cyan]{i}[/cyan]: {mod.name}")
    
    console.print(f"  [cyan]{len(modules) + 1}[/cyan]: Skip")

    choices = [str(i) for i in range(1, len(modules) + 2)]
    choice_str = Prompt.ask("Enter the number of your choice", choices=choices, default=str(len(modules) + 1))
    choice_idx = int(choice_str) - 1

    if choice_idx == len(modules):
        return None
    return modules[choice_idx]

run_app = typer.Typer(
    name="run",
    help="Run an agent system. Prompts for configuration if not provided via flags.",
    no_args_is_help=True,
)

class AppContext:
    def __init__(self):
        self.console = Console()
        self.agent_system: 'AgentSystem' | None = None
        self.driver_agent_name: str | None = None
        self.roster_instructions: str | None = None
        self.analysis_context: str | None = None
        self.sandbox_manager: 'SandboxManager' | None = None
        self.llm_client: object | None = None
        self.initial_history: List[dict] | None = None
        self.dataset_path: Path | None = None
        self.reference_dataset_path: Optional[Path] = None
        self.resources: List[Tuple[Path, str]] = []
        self.sandbox_details: dict = {}

@run_app.callback(invoke_without_command=True)
def main_run_callback(
    ctx: typer.Context,
    blueprint: Path = typer.Option(None, "--blueprint", "-bp", help="Path to the agent system JSON blueprint.", readable=True),
    driver_agent: str = typer.Option(None, "--driver-agent", "-d", help="Name of the agent to start with."),
    dataset: Path = typer.Option(None, "--dataset", "-ds", help="Path to the primary dataset file (.h5ad).", readable=True),
    reference_dataset: Path = typer.Option(None, "--reference-dataset", "-ref", help="Path to an optional reference dataset file (.h5ad).", readable=True),
    resources_dir: Path = typer.Option(None, "--resources", help="Path to a directory of resource files to mount.", exists=True, file_okay=False),
    llm_backend: str = typer.Option(None, "--llm", help="LLM backend to use: 'chatgpt' or 'ollama'."),
    ollama_host: str = typer.Option("http://localhost:11434", "--ollama-host", help="Base URL for Ollama backend."),
    sandbox: str = typer.Option(None, "--sandbox", help="Sandbox backend to use: 'docker' or 'singularity'."),
    force_refresh: bool = typer.Option(False, "--force-refresh", help="Force refresh/rebuild of the sandbox environment."),
):
    from olaf.agents.AgentSystem import AgentSystem
    from olaf.core.io_helpers import collect_resources
    from olaf.core.sandbox_management import init_docker, init_singularity_exec
    from olaf.datasets.czi_datasets import get_datasets_dir

    load_dotenv(dotenv_path=ENV_FILE)
    app_context = AppContext()
    console = app_context.console
    ctx.obj = app_context

    if blueprint is None:
        blueprint = _prompt_for_file(console, DEFAULT_AGENT_DIR, PACKAGE_AGENTS_DIR, ".json", "Agent System Blueprint")
    app_context.agent_system = AgentSystem.load_from_json(str(blueprint))
    
    if driver_agent is None:
        driver_agent = _prompt_for_driver(console, app_context.agent_system)
    if driver_agent not in app_context.agent_system.agents:
        raise typer.BadParameter(f"Driver agent '{driver_agent}' not found.")
    app_context.driver_agent_name = driver_agent
    app_context.roster_instructions = app_context.agent_system.get_instructions()
    
    if dataset is None:
        dataset = _prompt_for_file(console, get_datasets_dir(), PACKAGE_DATASETS_DIR, ".h5ad", "Primary Dataset")
    app_context.dataset_path = dataset

    if reference_dataset is None:
        if Prompt.ask("Do you want to add a reference dataset?", choices=["y", "n"], default="n").lower() == 'y':
            reference_dataset = _prompt_for_file(console, get_datasets_dir(), PACKAGE_DATASETS_DIR, ".h5ad", "Reference Dataset")
    app_context.reference_dataset_path = reference_dataset

    if sandbox is None:
        sandbox = Prompt.ask("Choose a sandbox backend", choices=["docker", "singularity"], default="docker")
    
    console.print(f"[cyan]Initializing sandbox backend: {sandbox}[/cyan]")
    script_dir = Path(__file__).resolve().parent
    
    manager_class, handle, copy_cmd, exec_endpoint, status_endpoint = (None, None, None, None, None)
    if sandbox == "docker":
        manager_class, handle, copy_cmd, exec_endpoint, status_endpoint = init_docker(script_dir, subprocess, console, force_refresh=force_refresh)
    elif sandbox == "singularity":
        manager_class, handle, copy_cmd, exec_endpoint, status_endpoint = init_singularity_exec(script_dir, SANDBOX_DATA_PATH, subprocess, console, force_refresh=force_refresh)
    else:
        raise typer.BadParameter(f"Unknown sandbox type '{sandbox}'. Supported: 'docker', 'singularity'.")
    app_context.sandbox_manager = manager_class()
    app_context.sandbox_details = {"handle": handle, "copy_cmd": copy_cmd, "is_exec_mode": sandbox == "singularity"}

    if llm_backend is None:
        llm_backend = Prompt.ask("Choose an LLM backend", choices=["chatgpt", "ollama"], default="chatgpt")
    if llm_backend == "ollama" and ollama_host == "http://localhost:11434":
         ollama_host = Prompt.ask("Enter the Ollama base URL", default="http://localhost:11434")

    console.print(f"[cyan]Initializing LLM backend: {llm_backend}[/cyan]")
    if llm_backend == "chatgpt":
        from openai import OpenAI
        app_context.llm_client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    elif llm_backend == "ollama":
        from olaf.core.ollama_wrapper import OllamaClient as OpenAI
        app_context.llm_client = OpenAI(host=ollama_host)
    else:
        raise typer.BadParameter(f"Unknown LLM backend '{llm_backend}'.")

    app_context.resources = collect_resources(console, resources_dir) if resources_dir else []
    if app_context.reference_dataset_path:
        app_context.resources.append((app_context.reference_dataset_path, SANDBOX_REF_DATA_PATH))

    analysis_context_str = f"Primary dataset path: **{SANDBOX_DATA_PATH}**\n"
    if app_context.reference_dataset_path:
        analysis_context_str += f"Reference dataset path: **{SANDBOX_REF_DATA_PATH}**\n"
    app_context.analysis_context = textwrap.dedent(analysis_context_str)
    
    driver = app_context.agent_system.get_agent(driver_agent)
    system_prompt = (app_context.roster_instructions + "\n\n" + driver.get_full_prompt(app_context.agent_system.global_policy) + "\n\n" + app_context.analysis_context)
    app_context.initial_history = [{"role": "system", "content": system_prompt}]

def _setup_and_run_session(context: AppContext, history: list, is_auto: bool, max_turns: int, benchmark_modules: Optional[List[Path]] = None):
    """Helper to start, run, and stop the sandbox session."""
    # --- Heavy imports needed for the session are deferred to here ---
    from olaf.execution.runner import run_agent_session, SandboxManager
    from olaf.agents.AgentSystem import AgentSystem
    from olaf.core.io_helpers import save_chat_history_as_json, save_chat_history_as_notebook

    sandbox_manager = cast(SandboxManager, context.sandbox_manager)
    console = context.console
    
    console.print("[cyan]Starting sandbox...[/cyan]")
    
    details = context.sandbox_details
    dataset_path = cast(Path, context.dataset_path)
    if details["is_exec_mode"] and hasattr(sandbox_manager, "set_data"):
        all_resources = [(dataset_path, SANDBOX_DATA_PATH)] + context.resources
        sandbox_manager.set_data(all_resources)
    if not sandbox_manager.start_container():
        console.print("[bold red]Failed to start sandbox container.[/bold red]")
       
```
