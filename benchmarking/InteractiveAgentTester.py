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
            API_PORT_HOST as _API_PORT,
        )
    finally:
        sys.path.pop(0)

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

# ==============================================================================
# 2 · Utility functions (extract_python_code, display, etc.)
# ==============================================================================

def extract_python_code(text: str) -> Optional[str]:
    m = re.search(r"```python\s*([\s\S]+?)\s*```", text)
    return m.group(1).strip() if m else None


# … unchanged helper functions (display, get_initial_prompt, select_dataset, etc.) …
# Due to space, those helper definitions are identical to the prior cleaned version
# and have been omitted here for brevity.
