#!/usr/bin/env python3
"""
Interactive Agent Tester – Docker, Singularity‑API, or **Singularity‑Exec (offline‑REPL)**
=======================================================================================
Run a natural‑language chat loop that generates runnable Python, executes it inside a
container, and streams the results back. Works even on clusters where **no networking**
is allowed for Singularity by using a long‑lived REPL inside the container.

Back‑ends
---------
1. **docker**            – Docker daemon + container with FastAPI kernel.
2. **singularity**       – Singularity *instance* with FastAPI kernel.
3. **singularity-exec**  – Long‑lived `singularity exec` REPL that talks to
                           `/opt/offline_kernel.py --repl` (no TCP).
"""
from __future__ import annotations

import base64
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import textwrap
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── 3rd‑party deps ──────────────────────────────────────────────────────────
try:
    from dotenv import load_dotenv
    from openai import OpenAI, APIError
    import requests  # only needed for networked back‑ends
    from rich.console import Console
    from rich.panel import Panel
    from rich.prompt import Prompt
    from rich.syntax import Syntax
    from rich.table import Table
except ImportError as e:
    print(f"Missing dependency: {e}. Install required packages.", file=sys.stderr)
    sys.exit(1)

console = Console()
SCRIPT_DIR = Path(__file__).resolve().parent
DATASETS_DIR = SCRIPT_DIR / "datasets"
OUTPUTS_DIR = SCRIPT_DIR / "outputs"
ENV_FILE = SCRIPT_DIR / ".env"

# In‑container canonical paths
SANDBOX_DATA_PATH = "/workspace/dataset.h5ad"
SANDBOX_RESOURCES_DIR = "/workspace/resources"

# ==============================================================================
# 1 · Choose back‑end BEFORE importing heavy managers
# ==============================================================================
backend = Prompt.ask(
    "Choose sandbox backend",
    choices=["docker", "singularity", "singularity-exec"],
    default="docker",
)

# Ask user whether to force‑update the sandbox image/SIF
force_refresh = (
    Prompt.ask(
        "Force update sandbox environment?", choices=["y", "n"], default="n"
    ).lower()
    == "y"
)

is_exec_mode = backend == "singularity-exec"

# -----------------------------------------------------------------------------
# 1a · Docker (FastAPI) back‑end
# -----------------------------------------------------------------------------
if backend == "docker":
    sandbox_dir = SCRIPT_DIR / "sandbox"
    sys.path.insert(0, str(sandbox_dir))
    try:
        from benchmarking_sandbox_management import (
            SandboxManager as _BackendManager,
            CONTAINER_NAME as _SANDBOX_HANDLE,
            IMAGE_NAME as _SANDBOX_IMAGE,  # assume this constant exists
            API_PORT_HOST as _API_PORT,
        )
    finally:
        sys.path.pop(0)

    # --- optional force‑refresh logic --------------------------------------
    if force_refresh:
        console.print("[yellow]Forcing Docker sandbox refresh…[/yellow]")
        # Stop & remove any running container gracefully
        subprocess.run(["docker", "rm", "-f", _SANDBOX_HANDLE], check=False)
        # Remove the sandbox image to ensure re‑pull/build
        subprocess.run(["docker", "image", "rm", "-f", _SANDBOX_IMAGE], check=False)
        console.print("[green]Docker image removed – it will be pulled/built on next start.[/green]")

    def COPY_CMD(src: str, dst: str):
        subprocess.run(["docker", "cp", src, dst], check=True)

    EXECUTE_ENDPOINT = f"http://localhost:{_API_PORT}/execute"
    STATUS_ENDPOINT = f"http://localhost:{_API_PORT}/status"

# -----------------------------------------------------------------------------
# 1b · Singularity instance (FastAPI) back‑end
# -----------------------------------------------------------------------------
elif backend == "singularity":
    sandbox_dir = SCRIPT_DIR / "sandbox"
    sys.path.insert(0, str(sandbox_dir))
    try:
        import benchmarking_sandbox_management_singularity as sing
    finally:
        sys.path.pop(0)

    # optional force‑refresh
    if force_refresh:
        console.print("[yellow]Forcing Singularity sandbox refresh…[/yellow]")
        try:
            sing.stop_instance()
        except Exception:
            pass  # ignore if not running
        if sing.SIF_PATH.exists():
            sing.SIF_PATH.unlink()
            console.print(
                f"[green]Deleted {sing.SIF_PATH.name} – it will be re‑downloaded on next start.[/green]"
            )

    class _SingInstanceWrapper:
        def start_container(self):
            return sing.start_instance()

        def stop_container(self):
            return sing.stop_instance()

    _BackendManager = _SingInstanceWrapper
    _SANDBOX_HANDLE = sing.INSTANCE_NAME
    _API_PORT = sing.API_PORT_HOST

    def COPY_CMD(src: str, dst: str):
        console.print(
            f"[yellow]Singularity instance: ensure {src} is reachable at {dst} via bind mount.[/yellow]"
        )

    EXECUTE_ENDPOINT = f"http://localhost:{_API_PORT}/execute"
    STATUS_ENDPOINT = f"http://localhost:{_API_PORT}/status"

# -----------------------------------------------------------------------------
# 1c · Singularity exec (offline‑REPL) back‑end
# -----------------------------------------------------------------------------
elif backend == "singularity-exec":
    sandbox_dir = SCRIPT_DIR / "sandbox"
    sys.path.insert(0, str(sandbox_dir))
    try:
        import benchmarking_sandbox_management_singularity as sing
    finally:
        sys.path.pop(0)

    # optional force‑refresh
    if force_refresh:
        console.print("[yellow]Forcing Singularity sandbox refresh…[/yellow]")
        if sing.SIF_PATH.exists():
            sing.SIF_PATH.unlink()
            console.print(
                f"[green]Deleted {sing.SIF_PATH.name} – it will be re‑downloaded on next start.[/green]"
            )

    SIF_PATH = sing.SIF_PATH
    SING_BIN = sing.SING_BIN
    SENTINEL = "<<<EOF>>>"

    class _SingExecBackend:
        """Launch one long‑lived REPL inside the SIF and stream code to it."""

        def __init__(self):
            self._binds: List[str] = []
            self._proc: Optional[subprocess.Popen[str]] = None

        def set_data(self, dataset: Path, resources: List[Tuple[Path, str]]):
            self._binds = [
                "--bind",
                f"{dataset.resolve()}:{SANDBOX_DATA_PATH}",
            ]
            for host, cont in resources:
                self._binds.extend(["--bind", f"{host.resolve()}:{cont}"])

        # ------------------------------------------------------------------
        # Container lifecycle
        # ------------------------------------------------------------------
        def start_container(self):
            if self._proc:
                return True  # already running
            if not sing.pull_sif_if_needed():
                return False

            cmd = [
                SING_BIN,
                "exec",
                "--containall",
                "--cleanenv",
                *self._binds,
                str(SIF_PATH),
                "python",
                "/opt/offline_kernel.py",
                "--repl",
            ]
            self._proc = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,  # line buffered
            )
            # Wait for the REPL banner
            ready_line = self._proc.stdout.readline().strip()
            if ready_line != "__REPL_READY__":
                console.print(
                    f"[red]REPL failed to start. Got: {ready_line}[/red]"
                )
                self.stop_container()
                return False
            return True

        def stop_container(self):
            if not self._proc:
                return True
            try:
                if self._proc.stdin:
                    self._proc.stdin.close()
                self._proc.terminate()
                self._proc.wait(timeout=5)
            except Exception:
                self._proc.kill()
            self._proc = None
            return True

        # ------------------------------------------------------------------
        # Code execution
        # ------------------------------------------------------------------
        def exec_code(self, code: str, timeout: int = 300) -> Dict:
            if not self._proc:
                raise RuntimeError("REPL not running")
            assert self._proc.stdin and self._proc.stdout

            # Send code block + sentinel
            self._proc.stdin.write(code)
            if not code.endswith("\n"):
                self._proc.stdin.write("\n")
            self._proc.stdin.write(SENTINEL + "\n")
            self._proc.stdin.flush()

            # Read exactly one JSON line
            start_time = time.time()
            while True:
                if time.time() - start_time > timeout:
                    return {
                        "status": "timeout",
                        "stdout": "",
                        "stderr": "Execution timed out in REPL.",
                        "images": [],
                    }
                line = self._proc.stdout.readline()
                if not line:
                    continue
                line = line.strip()
                try:
                    return json.loads(line)
                except json.JSONDecodeError:
                    # Non‑JSON noise; continue reading
                    continue

    _BackendManager = _SingExecBackend

    def COPY_CMD(src: str, dst: str):
        console.print("[yellow]singularity-exec mode uses bind mounts instead of docker cp.[/yellow]")
else:
    console.print("[red]Unknown backend.")
    sys.exit(1)

# ====================================================================================
# 2 · Generic helpers (unchanged)
# ====================================================================================

def extract_python_code(txt: str) -> Optional[str]:
    m = re.search(r"```python\s*([\s\S]+?)\s*```", txt)
    return m.group(1).strip() if m else None


# Rich display wrappers

def _panel(role: str, content: str):
    titles = {"system": "SYSTEM", "user": "USER", "assistant": "ASSISTANT"}
    styles = {"system": "dim blue", "user": "cyan", "assistant": "green"}
    console.print(Panel(content, title=titles.get(role, role.upper()), border_style=styles.get(role, "white")))


def display(role: str, content: str):
    if role == "assistant":
        code = extract_python_code(content) or ""
        text_part = re.sub(r"```python[\s\S]+?```", "", content, count=1).strip()
        if text_part:
            _panel("assistant", text_part)
        if code:
            console.print(
                Panel(
                    Syntax(code, "python", line_numbers=True),
                    title="ASSISTANT (code)",
                    border_style="green",
                )
            )
    else:
        _panel(role, content)


# ====================================================================================
# 3 · Dataset / prompt helpers (unchanged)
# ====================================================================================

def get_initial_prompt() -> str:
    console.print("[bold cyan]Enter the initial user prompt (Ctrl+D to finish):[/bold cyan]")
    try:
        txt = sys.stdin.read().strip()
    except EOFError:
        txt = ""
    if not txt:
        console.print("[red]Empty prompt – aborting.[/red]")
        sys.exit(1)
    return txt


def select_dataset() -> Tuple[Path, dict]:
    if not DATASETS_DIR.exists():
        console.print(f"[red]Datasets dir not found: {DATASETS_DIR}[/red]")
        sys.exit(1)
    items = [
        (p, json.loads(p.with_suffix(".json").read_text()))
        for p in DATASETS_DIR.glob("*.h5ad")
        if p.with_suffix(".json").exists()
    ]
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
    idx = int(Prompt.ask("Choose index", choices=[str(i) for i in range(1, len(items) + 1)])) - 1
    return items[idx]


def collect_resources() -> List[Tuple[Path, str]]:
    console.print("\n[bold cyan]Optional: paths to bind inside sandbox[/bold cyan] (blank line to finish)")
    res: List[Tuple[Path, str]] = []
    while True:
        p = Prompt.ask("Path", default="").strip()
        if not p:
            break
        path = Path(p).expanduser().resolve()
        if not path.exists():
            console.print(f"[yellow]Path does not exist: {path}[/yellow]")
            continue
        res.append((path, f"{SANDBOX_RESOURCES_DIR}/{path.name}"))
    return res


# ====================================================================================
# 4 · Networked FastAPI helpers (skipped for exec mode)
# ====================================================================================

def api_alive(max_retries: int = 10, delay: float = 1.5) -> bool:
    if is_exec_mode:
        return True  # nothing to ping
    for _ in range(max_retries):
        try:
            if requests.get(STATUS_ENDPOINT, timeout=2).json().get("status") == "ok":
                return True
        except Exception:
            time.sleep(delay)
    return False


def format_execute_response(resp: dict) -> str:
    lines = ["Code execution result:"]
    if resp.get("status") != "ok":
        lines.append(f"[status: {resp.get('status')}]")
    stdout, stderr = resp.get("stdout", ""), resp.get("stderr", "")
    if stdout:
        lines += ["--- STDOUT ---", stdout[:1500]]
    if stderr:
        lines += ["--- STDERR ---", stderr[:1500]]
    img_paths = []
    for b64 in resp.get("images", []):
        fname = OUTPUTS_DIR / f"{datetime.now():%Y%m%d_%H%M%S_%f}.png"
        fname.parent.mkdir(exist_ok=True, parents=True)
        with open(fname, "wb") as f:
            f.write(base64.b64decode(b64))
        img_paths.append(str(fname))
    if img_paths:
        lines.append("Saved images: " + ", ".join(img_paths))
    return "\n".join(lines)


# ====================================================================================
# 5 · Main interactive loop (unchanged)
# ====================================================================================

def run_interactive(prompt: str, dataset: Path, metadata: dict, resources: List[Tuple[Path, str]]):
    mgr = _BackendManager()
    console.print(f"Starting sandbox ({backend}) …")

    # Tell exec back‑end where data/resources are (creates bind list)
    if is_exec_mode and hasattr(mgr, "set_data"):
        mgr.set_data(dataset, resources)

    if not mgr.start_container():
        console.print("[red]Failed to start sandbox.[/red]")
        return

    if not api_alive():
        console.print("[red]Kernel API not responsive (networked back‑end).[/red]")
        return

    # For docker / singularity‑instance we still *attempt* docker cp (no‑op or warning otherwise)
    if not is_exec_mode:
        COPY_CMD(str(dataset), f"{_SANDBOX_HANDLE}:{SANDBOX_DATA_PATH}")
        for h, c in resources:
            COPY_CMD(str(h), f"{_SANDBOX_HANDLE}:{c}")

    resource_lines = [f"- {c} (from {h})" for h, c in resources] or ["- (none)"]
    sys_prompt = textwrap.dedent(
        f"""
        You are an AI assistant analysing a single‑cell dataset.
        Dataset path inside container: **{SANDBOX_DATA_PATH}**
        Additional resources:\n"""
        + "\n".join(resource_lines)
        + "\n\n"
        + textwrap.dedent(
            f"Dataset metadata:\n{json.dumps(metadata, indent=2)}\n\n"
            "Wrap runnable Python in triple‑backtick ```python blocks. Imports & variables persist within the container session."
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
            rsp = openai.chat.completions.create(
                model="gpt-4o", messages=history, temperature=0.7
            )
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
                if is_exec_mode:
                    exec_result = mgr.exec_code(code, timeout=300)
                else:
                    exec_result = requests.post(
                        EXECUTE_ENDPOINT, json={"code": code, "timeout": 300}, timeout=310
                    ).json()
                feedback = format_execute_response(exec_result)
            except Exception as exc:
                feedback = f"Code execution result:\n[Execution error on host: {exc}]"

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

    console.print("Stopping sandbox…")
    mgr.stop_container()


# ====================================================================================
# 6 · Entry‑point
# ====================================================================================

def main():
    load_dotenv(Path(ENV_FILE))
    if not os.getenv("OPENAI_API_KEY"):
        console.print(f"[red]OPENAI_API_KEY not set in {ENV_FILE}.[/red]")
        sys.exit(1)

    prompt = get_initial_prompt()
    data_p, meta = select_dataset()
    resources = collect_resources()
    run_interactive(prompt, data_p, meta, resources)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        console.print("\nInterrupted.")
