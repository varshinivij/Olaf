#!/usr/bin/env python3
"""
Run arbitrary Python and emit a JSON blob describing the result.
Designed for `singularity exec` with no network.
"""
from __future__ import annotations
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/.matplotlib")

import base64
import io
import json
import sys
import traceback
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

# If you want stateful sessions later, replace this with a long‐lived namespace.
global_namespace: dict = {}

def run_code(code: str, namespace: dict | None = None) -> dict:
    """
    Execute a Python code string, capturing stdout, stderr, and any Matplotlib figures.

    Args:
        code: Python source to execute.
        namespace: Optional dict for globals() during execution. If None, a new dict is used.

    Returns:
        A dict with keys "status", "stdout", "stderr", and "images" (base64 PNGs).
    """
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    images: list[str] = []
    status = "ok"

    ns = namespace if namespace is not None else {}
    ns.setdefault("__builtins__", __builtins__)

    try:
        with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
            # Ensure Matplotlib uses Agg backend if available
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                ns["plt"] = plt
            except ImportError:
                pass

            compiled = compile(code, "<user_code>", "exec")
            exec(compiled, ns)

            # Capture open figures
            if "plt" in ns:
                for fig_num in ns["plt"].get_fignums():
                    fig = ns["plt"].figure(fig_num)
                    buf = io.BytesIO()
                    fig.savefig(buf, format="png", bbox_inches="tight")
                    images.append(base64.b64encode(buf.getvalue()).decode())
                    ns["plt"].close(fig)
    except Exception:
        stderr_buf.write(traceback.format_exc())
        status = "error"
    finally:
        if "plt" in ns:
            try:
                ns["plt"].close("all")
            except Exception:
                pass

    return {
        "status": status,
        "stdout": stdout_buf.getvalue(),
        "stderr": stderr_buf.getvalue(),
        "images": images,
    }

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: offline_kernel.py <code_file_path>", file=sys.stderr)
        sys.exit(2)

    code_path = Path(sys.argv[1])
    if not code_path.is_file():
        result = {
            "status": "error",
            "stdout": "",
            "stderr": f"Error: Code file '{code_path}' not found.",
            "images": [],
        }
        print(json.dumps(result))
        sys.exit(3)

    try:
        user_code = code_path.read_text(encoding="utf-8")
    except Exception as e:
        result = {
            "status": "error",
            "stdout": "",
            "stderr": f"Error reading '{code_path}': {e}",
            "images": [],
        }
        print(json.dumps(result))
        sys.exit(4)

    # Each invocation uses the same global_namespace (stateless for now).
    output = run_code(user_code, namespace=global_namespace)
    print(json.dumps(output))