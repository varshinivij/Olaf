#!/usr/bin/env python3
"""
Interactive Agent Tester (Docker **or** Singularity backend)
==========================================================
A unified interactive tester that can drive either the **Docker sandbox** (`benchmarking_sandbox_management.py`)
or the **Apptainer/Singularity sandbox** (`benchmarking_sandbox_management_singularity.py`).

At launch you choose a backend:
    • *docker*       – requires Docker daemon on this machine.
    • *singularity*  – requires `apptainer`/`singularity`; no Docker needed.

The rest of the behaviour (multi‑turn GPT orchestration, FastAPI kernel execution,
resource upload, unlimited chat loop) is unchanged.
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

# ── Third‑party deps ─────────────────────────────────────────────────────────
try:
    from dotenv import load_dotenv
    from openai import OpenAI, APIError
    import requests
    from rich.console import Console
    from rich.panel import Panel
    from rich.prompt import Prompt
    from rich.syntax import Syntax
    from rich.table import Table
except ImportError as e:
    print(f"Missing dependency: {e}. Install required packages.", file=sys.stderr)
    sys.exit(1)

console = Console()

# ── Runtime‑backend selection (ask the user **before** importing managers) ──
backend = Prompt.ask("Choose sandbox backend", choices=["docker", "singularity"], default="docker")

SCRIPT_DIR = Path(__file__).resolve().parent

if backend == "docker":
    sandbox_dir = SCRIPT_DIR / "sandbox"
    sys.path.insert(0, str(sandbox_dir))
    try:
        from benchmarking_sandbox_management import (
            SandboxManager as _BackendManager,
            CONTAINER_NAME as _SANDBOX_HANDLE,
            API_PORT_HOST as _API_PORT,
        )
    finally:
        sys.path.pop(0)
    COPY_CMD = lambda src, dst: subprocess.run(["docker", "cp", src, dst], check=True)

elif backend == "singularity":
    sandbox_dir = SCRIPT_DIR / "sandbox"
    sys.path.insert(0, str(sandbox_dir))
    try:
        import benchmarking_sandbox_management_singularity as sing
    except ImportError as e:
        console.print(f"[red]Failed to import Singularity manager: {e}[/red]")
        sys.exit(1)

    class _SingWrapper:  # thin adapter to mimic Docker SandboxManager API
        def __init__(self):
            pass
        def start_container(self):
            return sing.start_instance()
        def stop_container(self, remove: bool = True, container_obj=None):
            return sing.stop_instance()
    _BackendManager = _SingWrapper
    _SANDBOX_HANDLE = sing.INSTANCE_NAME
    _API_PORT = sing.API_PORT_HOST

    # Apptainer/ Singularity lacks a simple cp, so we issue a warning and rely on bind‑mounts
    def COPY_CMD(src, dst):  # noqa: N802
        console.print(f"[yellow]File copy inside Singularity instance not automated.\n"
                      f"Ensure the file {src} is reachable at {dst} via bind mount or in the definition file.[/yellow]")

else:
    console.print("[red]Unknown backend choice.[/red]")
    sys.exit(1)

# ── Constants (after backend choice) ─────────────────────────────────────────
DATASETS_DIR = SCRIPT_DIR / "datasets"
OUTPUTS_DIR = SCRIPT_DIR / "outputs"
ENV_FILE = SCRIPT_DIR / ".env"
SANDBOX_DATA_PATH = "/home/sandboxuser/data.h5ad"
SANDBOX_RESOURCES_DIR = "/home/sandboxuser/resources"
API_BASE_URL = f"http://localhost:{_API_PORT}"
EXECUTE_ENDPOINT = f"{API_BASE_URL}/execute"
STATUS_ENDPOINT = f"{API_BASE_URL}/status"


# ── Helper utilities ────────────────────────────────────────────────────────

def extract_python_code(txt: str) -> str | None:
    m = re.search(r"```python\s*([\s\S]+?)\s*```", txt)
    return m.group(1).strip() if m else None


def display(role: str, content: str) -> None:
    titles = {"system": "SYSTEM", "user": "USER", "assistant": "ASSISTANT"}
    styles = {"system": "dim blue", "user": "cyan", "assistant": "green"}
    title = titles.get(role, role.upper())
    style = styles.get(role, "white")

    if role == "assistant":
        code = extract_python_code(content)
        txt = re.sub(r"```python[\s\S]+?```", "", content, count=1).strip()
        if txt:
            console.print(Panel(txt, title=f"{title} (text)", border_style=style))
        if code:
            console.print(Panel(Syntax(code, "python", line_numbers=True), title=f"{title} (code)", border_style=style))
    else:
        console.print(Panel(content, title=title, border_style=style))


# ── Dataset & prompt helpers ────────────────────────────────────────────────

def get_initial_prompt() -> str:
    console.print("[bold cyan]Enter the initial user prompt (Ctrl+D to finish):[/bold cyan]")
    try:
        txt = sys.stdin.read().strip()
    except EOFError:
        txt = ""
    if not txt:
        console.print("[red]Empty prompt. Aborting.[/red]")
        sys.exit(1)
    return txt


def select_dataset() -> Tuple[Path, dict]:
    if not DATASETS_DIR.exists():
        console.print(f"[red]Datasets dir not found: {DATASETS_DIR}[/red]")
        sys.exit(1)
    items = [(p, json.loads(p.with_suffix(".json").read_text())) for p in DATASETS_DIR.glob("*.h5ad") if p.with_suffix(".json").exists()]
    if not items:
        console.print("[red]No datasets found.[/red]")
        sys.exit(1)
    tbl = Table(title="Datasets")
    tbl.add_column("Idx", justify="right")
    tbl.add_column("Name")
    tbl.add_column("Cells", justify="right")
    for i, (p, meta) in enumerate(items, 1):
        tbl.add_row(str(i), meta.get("dataset_title", p.stem), str(meta.get("cell_count", "?")))
    console.print(tbl)
    idx = int(Prompt.ask("Choose index", choices=[str(i) for i in range(1, len(items)+1)])) - 1
    return items[idx]


def collect_resources() -> List[Tuple[Path, str]]:
    console.print("\n[bold cyan]Optional: list files/folders to copy into sandbox[/bold cyan] (blank line to finish)")
    lst: List[Tuple[Path, str]] = []
    while True:
        p = Prompt.ask("Path", default="").strip()
        if not p:
            break
        path = Path(p).expanduser().resolve()
        if not path.exists():
            console.print(f"[yellow]Path does not exist: {path}[/yellow]")
            continue
        lst.append((path, f"{SANDBOX_RESOURCES_DIR}/{path.name}"))
    return lst


# ── FastAPI kernel helpers ──────────────────────────────────────────────────

def api_alive(max_retries: int = 10, delay: float = 1.5) -> bool:
    for _ in range(max_retries):
        try:
            if requests.get(STATUS_ENDPOINT, timeout=2).json().get("status") == "ok":
                return True
        except requests.RequestException:
            time.sleep(delay)
    return False


def format_execute_response(resp: dict) -> str:
    lines = ["Code execution result:"]
    stdout, stderr, imgs = [], [], []
    for itm in resp.get("outputs", []):
        if itm["type"] == "stream":
            (stdout if itm.get("name") == "stdout" else stderr).append(itm.get("text", ""))
        elif itm["type"] == "error":
            stderr.append("Error: " + itm.get("evalue", ""))
            stderr.extend(itm.get("traceback", []))
        elif itm["type"] == "display_data":
            for mime, b64 in itm.get("data", {}).items():
                if mime.startswith("image/"):
                    fname = OUTPUTS_DIR / f"{datetime.now():%Y%m%d_%H%M%S_%f}.{mime.split('/')[1].split('+')[0]}"
                    fname.parent.mkdir(exist_ok=True)
                    with open(fname, "wb") as f:
                        f.write(base64.b64decode(b64))
                    imgs.append(str(fname))
    if stdout:
        lines += ["--- STDOUT ---", "".join(stdout)[:1500]]
    if stderr:
        lines += ["--- STDERR ---", "".join(stderr)[:1500]]
    if imgs:
        lines.append("Saved images: " + ", ".join(imgs))
    lines.append(f"Final Status: {resp.get('final_status')}")
    return "\n".join(lines)


# ── Chat‑runner ─────────────────────────────────────────────────────────────

def run_interactive(prompt: str, dataset: Path, metadata: dict, resources: List[Tuple[Path, str]]) -> None:
    mgr = _BackendManager()
    console.print(f"Starting sandbox ({backend}) …")
    if not mgr.start_container():
        console.print("[red]Failed to start sandbox.[/red]")
        return

    try:
        if not api_alive():
            console.print("[red]Kernel API not responsive.[/red]")
            return
        # dataset copy (Docker only, Singularity warns via COPY_CMD)
        COPY_CMD(str(dataset), f"{_SANDBOX_HANDLE}:{SANDBOX_DATA_PATH}")
        for h, c in resources:
            COPY_CMD(str(h), f"{_SANDBOX_HANDLE}:{c}")

        resource_lines = [f"- {c} (from {h})" for h, c in resources] or ["- (none)"]
        sys_prompt = textwrap.dedent(
            f"""
            You are an AI assistant analysing a single‑cell dataset.  The file lives inside the sandbox at **{SANDBOX_DATA_PATH}**.
            Additional resources:\n""" + "\n".join(resource_lines) + "\n\n" + textwrap.dedent(
                f"Dataset metadata:\n{json.dumps(metadata, indent=2)}\n\nWrap runnable Python in triple‑backtick ```python blocks. Imports & vars persist."""
            )
        )

        history = [
            {"role": "system", "content": sys_prompt},
            {"role": "user", "content": prompt},
        ]
        display("system", sys_prompt)
        display("user", prompt)

        openai = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))
        turn = 0
        while True:
            turn += 1
            console.print(f"\n[bold]OpenAI call (turn {turn})…[/bold]")
            try:
                rsp = openai.chat.completions.create(model="gpt-4o", messages=history, temperature=0.7)
            except APIError as e:
                console.print(f"[red]OpenAI error: {e}[/red]")
                break
            assistant_msg = rsp.choices[0].message.content
            history.append({"role": "assistant", "content": assistant_msg})
            display("assistant", assistant_msg)

            code = extract_python_code(assistant_msg)
            if code:
                console.print("[cyan]Executing code…[/cyan]")
                try:
                    api_r = requests.post(EXECUTE_ENDPOINT, json={"code": code, "timeout": 120}, timeout=130).json()
                    feedback = format_execute_response(api_r)
                except Exception as exc:
                    feedback = f"Code execution result:\n[Execution error: {exc}]"
                history.append({"role": "user", "content": feedback})
                display("user", feedback)

            console.print("\n[bold]Next message (blank = continue, 'exit' to quit):[/bold]")
            try:
                user_in = input().strip()
            except (EOFError, KeyboardInterrupt):
                user_in = "exit"
            if user_in.lower() in {"exit", "quit"}:
                break
            if user_in:
                history.append({"role": "user", "content": user_in})
                display("user", user_in)
    finally:
        console.print("Stopping sandbox…")
        mgr.stop_container(remove=True)


# ── CLI entry ───────────────────────────────────────────────────────────────

def main():
    load_dotenv(Path(ENV_FILE))
    if not os.getenv("OPENAI_API_KEY"):
        console.print(f"[red]OPENAI_API_KEY not set in {ENV_FILE}.[/red]")
        sys.exit(1)

    prompt = get_initial_prompt()
    data_p, meta = select_dataset()
    res = collect_resources()
    run_interactive(prompt, data_p, meta, res)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        console.print("\nInterrupted.")
