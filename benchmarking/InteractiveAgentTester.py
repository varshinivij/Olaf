#!/usr/bin/env python3
"""
Interactive Agent Tester (API Mode)
==================================
A drop‑in replacement for the previous *one‑shot* tester.  This version keeps the
same execution model (GPT‑powered assistant that sends code which is executed in
an isolated Docker sandbox exposing a FastAPI kernel service) but removes the
hard limit on message turns and lets **you** steer the dialogue interactively.

Key additions
-------------
* **Unlimited conversation** – after every turn you can type a follow‑up message
  that is appended to the chat history before the next assistant call.
* **Resource upload** – copy any number of files or whole folders from the host
  into the running sandbox (stored under */home/sandboxuser/resources/*).  A
  summary of the uploaded paths is automatically prepended to the system prompt
  so the assistant knows what is available.

Usage
-----
$ python interactive_agent_tester.py  # guided TUI (Rich)

Exit the interactive loop at any time by typing `exit`, `quit` or pressing
`Ctrl‑C`.
"""

from __future__ import annotations

import argparse
import base64
import json
import os
import re
import shlex
import subprocess
import sys
import textwrap
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

# --- Third‑party deps ---------------------------------------------------------
try:
    from dotenv import load_dotenv
    from openai import OpenAI, APIError
    import requests
    from rich.console import Console
    from rich.markdown import Markdown
    from rich.panel import Panel
    from rich.prompt import Prompt, Confirm
    from rich.syntax import Syntax
    from rich.table import Table
except ImportError as e:  # graceful fallback if Rich not installed
    print(f"Missing dependency: {e}.  Please install required packages.", file=sys.stderr)
    sys.exit(1)

# --- Local sandbox manager ----------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
sandbox_dir = SCRIPT_DIR / "sandbox"
sys.path.insert(0, str(sandbox_dir))
try:
    from benchmarking_sandbox_management import (
        SandboxManager,
        CONTAINER_NAME as SANDBOX_CONTAINER_NAME,
        API_PORT_HOST,
    )
finally:
    sys.path.pop(0)

# --- Constants ----------------------------------------------------------------
DATASETS_DIR = SCRIPT_DIR / "datasets"
OUTPUTS_DIR = SCRIPT_DIR / "outputs"
ENV_FILE = SCRIPT_DIR / ".env"
SANDBOX_DATA_PATH = "/home/sandboxuser/data.h5ad"
SANDBOX_RESOURCES_DIR = "/home/sandboxuser/resources"
API_BASE_URL = f"http://localhost:{API_PORT_HOST}"
EXECUTE_ENDPOINT = f"{API_BASE_URL}/execute"
STATUS_ENDPOINT = f"{API_BASE_URL}/status"
console = Console()

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

def extract_python_code(text: str) -> str | None:
    """Return the first ```python``` code block found in *text*, or *None*."""
    m = re.search(r"```python\s*([\s\S]+?)\s*```", text)
    return m.group(1).strip() if m else None


def display(role: str, content: str) -> None:
    """Pretty print a chat turn to the terminal using Rich formatting."""
    title_map = {
        "system": "SYSTEM",
        "user": "USER",
        "assistant": "ASSISTANT",
    }
    style_map = {
        "system": "dim blue",
        "user": "cyan",
        "assistant": "green",
    }
    title = title_map.get(role, role.upper())
    style = style_map.get(role, "white")

    if role == "assistant":
        code = extract_python_code(content)
        text_part = re.sub(r"```python[\s\S]+?```", "", content, count=1).strip()
        if text_part:
            console.print(Panel(text_part, title=f"{title} (text)", border_style=style))
        if code:
            console.print(Panel(Syntax(code, "python", line_numbers=True), title=f"{title} (code)", border_style=style))
    else:
        console.print(Panel(content, title=title, border_style=style))


# -----------------------------------------------------------------------------
# Prompts & selection helpers (unchanged except for small tweaks)
# -----------------------------------------------------------------------------

def get_initial_prompt() -> str:
    console.print("[bold cyan]Enter the initial user prompt for the agent.[/bold cyan]")
    console.print("Finish with Ctrl+D (Unix) / Ctrl+Z (Windows).")
    try:
        text = sys.stdin.read().strip()
    except EOFError:
        text = ""
    if not text:
        console.print("[red]Empty prompt. Aborting.[/red]")
        sys.exit(1)
    return text


def select_dataset() -> Tuple[Path, dict]:
    if not DATASETS_DIR.is_dir():
        console.print(f"[red]Datasets directory not found: {DATASETS_DIR}[/red]")
        sys.exit(1)
    datasets = []
    for p in DATASETS_DIR.glob("*.h5ad"):
        meta_path = p.with_suffix(".json")
        if meta_path.exists():
            datasets.append((p, json.loads(meta_path.read_text())))
    if not datasets:
        console.print("[red]No datasets found.[/red]")
        sys.exit(1)
    table = Table(title="Available datasets")
    table.add_column("Idx", justify="right")
    table.add_column("Name")
    table.add_column("Cells", justify="right")
    for i, (p, meta) in enumerate(datasets, 1):
        table.add_row(str(i), meta.get("dataset_title", p.stem), str(meta.get("cell_count", "?")))
    console.print(table)
    idx = int(Prompt.ask("Choose dataset index", choices=[str(i) for i in range(1, len(datasets) + 1)])) - 1
    return datasets[idx]


def collect_resources() -> List[Tuple[Path, str]]:
    """Prompt user for files/folders to add to sandbox.  Returns list of tuples
    (host_path, container_path)."""
    resources: List[Tuple[Path, str]] = []
    console.print("\n[bold cyan]Add extra resources to the sandbox (optional).[/bold cyan]")
    console.print("Enter absolute or relative paths one per line.  Leave empty line to finish.")
    while True:
        path_str = Prompt.ask("Path", default="").strip()
        if not path_str:
            break
        path = Path(path_str).expanduser().resolve()
        if not path.exists():
            console.print(f"[yellow]Path does not exist: {path}[/yellow]")
            continue
        container_dst = f"{SANDBOX_RESOURCES_DIR}/{path.name}"
        resources.append((path, container_dst))
    return resources


# -----------------------------------------------------------------------------
# API helpers
# -----------------------------------------------------------------------------

def api_alive(max_retries: int = 10, delay: float = 1.5) -> bool:
    for _ in range(max_retries):
        try:
            if requests.get(STATUS_ENDPOINT, timeout=2).json().get("status") == "ok":
                return True
        except requests.RequestException:
            time.sleep(delay)
    return False


def format_execute_response(resp: dict) -> str:
    out_lines = ["Code execution result:"]
    std_out, std_err = [], []
    images = []
    for item in resp.get("outputs", []):
        if item["type"] == "stream":
            (std_out if item.get("name") == "stdout" else std_err).append(item.get("text", ""))
        elif item["type"] == "error":
            std_err.append("Error: " + item.get("evalue", ""))
            std_err.extend(item.get("traceback", []))
        elif item["type"] == "display_data":
            for mime, b64 in item.get("data", {}).items():
                if mime.startswith("image/"):
                    fname = OUTPUTS_DIR / f"{datetime.now():%Y%m%d_%H%M%S_%f}.{mime.split('/')[1].split('+')[0]}"
                    fname.parent.mkdir(exist_ok=True)
                    with open(fname, "wb") as fh:
                        fh.write(base64.b64decode(b64))
                    images.append(str(fname))
    if std_out:
        out_lines.append("--- STDOUT ---")
        out_lines.append("".join(std_out)[:1500])
    if std_err:
        out_lines.append("--- STDERR ---")
        out_lines.append("".join(std_err)[:1500])
    if images:
        out_lines.append("Saved images: " + ", ".join(images))
    out_lines.append(f"Final Status: {resp.get('final_status')}")
    return "\n".join(out_lines)


# -----------------------------------------------------------------------------
# Main interactive runner
# -----------------------------------------------------------------------------

def run_interactive(prompt: str, dataset: Path, metadata: dict, resources: List[Tuple[Path, str]]) -> None:
    # 1. Start sandbox container
    mgr = SandboxManager()
    console.print("Starting sandbox container …")
    if not mgr.start_container():
        console.print("[red]Failed to start container.[/red]")
        return

    try:
        # 2. Wait for kernel API
        if not api_alive():
            console.print("[red]Kernel API did not become responsive.[/red]")
            return

        # 3. Copy dataset
        subprocess.run(["docker", "cp", str(dataset), f"{SANDBOX_CONTAINER_NAME}:{SANDBOX_DATA_PATH}"], check=True)

        # 4. Copy extra resources
        for host_path, cont_path in resources:
            subprocess.run(["docker", "cp", str(host_path), f"{SANDBOX_CONTAINER_NAME}:{cont_path}"], check=True)

        # 5. Build system prompt
        resource_lines = [f"- {cpath} (from {hpath})" for hpath, cpath in resources] or ["- (none)"]
        system_prompt = textwrap.dedent(
            f"""
            You are an AI assistant tasked with analysing a single‑cell transcriptomics dataset.
            The dataset is available at **{SANDBOX_DATA_PATH}** inside the execution environment.

            Additional resources copied for this session:\n""" + "\n".join(resource_lines) + "\n\n" + textwrap.dedent(
                f"""
                Dataset metadata:\n{json.dumps(metadata, indent=2)}

                Always wrap executable Python in a single triple‑backtick block with the language spec *python*.
                Variables and imports persist between executions.
            """
            )
        )

        # 6. Chat loop
        history = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": prompt},
        ]

        display("system", system_prompt)
        display("user", prompt)

        openai = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
        attempt = 0
        while True:
            attempt += 1
            console.print(f"\n[bold]OpenAI call (turn {attempt})…[/bold]")
            try:
                resp = openai.chat.completions.create(
                    model="gpt-4o", messages=history, temperature=0.7
                )
            except APIError as e:
                console.print(f"[red]OpenAI error: {e}[/red]")
                break

            assistant_msg = resp.choices[0].message.content
            history.append({"role": "assistant", "content": assistant_msg})
            display("assistant", assistant_msg)

            # Execute any code
            code = extract_python_code(assistant_msg)
            if code:
                console.print("[cyan]Executing code inside sandbox…[/cyan]")
                try:
                    api_resp = requests.post(EXECUTE_ENDPOINT, json={"code": code, "timeout": 120}, timeout=130).json()
                    user_feedback = format_execute_response(api_resp)
                except Exception as e:
                    user_feedback = f"Code execution result:\n[Execution error: {e}]"
                history.append({"role": "user", "content": user_feedback})
                display("user", user_feedback)

            # Ask user for next input
            console.print("\n[bold]Enter next message (blank to continue, 'exit' to quit):[/bold]")
            try:
                user_input = input().strip()
            except (EOFError, KeyboardInterrupt):
                user_input = "exit"
            if user_input.lower() in {"exit", "quit"}:
                console.print("[green]Ending session.[/green]")
                break
            if user_input:
                history.append({"role": "user", "content": user_input})
                display("user", user_input)
            # else: blank → assistant continues next loop

    finally:
        console.print("Stopping sandbox…")
        mgr.stop_container(remove=True)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> None:
    load_dotenv(dotenv_path=ENV_FILE)
    if not os.getenv("OPENAI_API_KEY"):
        console.print(f"[red]OPENAI_API_KEY not found in {ENV_FILE}.[/red]")
        sys.exit(1)

    console.print("[bold blue]Interactive Agent Tester (API Mode)[/bold blue]")
    prompt = get_initial_prompt()
    dataset_path, metadata = select_dataset()
    extra_resources = collect_resources()
    run_interactive(prompt, dataset_path, metadata, extra_resources)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        console.print("\nInterrupted. Goodbye.")
