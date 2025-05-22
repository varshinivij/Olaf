#!/usr/bin/env python3
"""Singularity Sandbox Manager (Docker‑free)
==========================================
Pure‑Singularity version that **does not require Docker at all**.  It expects a
`Singularity` (definition) file in the same directory and builds a `sandbox.sif`
from it when necessary.

Commands (same as before)
-------------------------
    build   – build `sandbox.sif` from the local Singularity file
    start   – start an *instance* exposing the FastAPI kernel on host port 8000
    stop    – stop & remove the instance
    status  – show instance + port status
    logs    – tail the instance log (default 50 lines)

Run with no args for an interactive REPL.
"""
from __future__ import annotations

import argparse
import logging
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DEF_FILE = SCRIPT_DIR / "Singularity"          # definition file expected here
SIF_PATH = SCRIPT_DIR / "sandbox.sif"          # output image
INSTANCE_NAME = "benchmarking_sandbox_instance"
API_PORT_INSIDE = 8000
API_PORT_HOST = 8000

SING_BIN = shutil.which("apptainer") or shutil.which("singularity")
if not SING_BIN:
    print("Singularity/Apptainer executable not found in PATH. Do you need to load a module?", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: list[str], *, capture: bool = False, check: bool = True):
    logging.debug("$ %s", " ".join(shlex.quote(c) for c in cmd))
    return subprocess.run(cmd, text=True, capture_output=capture, check=check)


def build_sif(rebuild: bool = False) -> bool:
    """Build sandbox.sif from local Singularity def file if needed."""
    if not DEF_FILE.exists():
        logging.error("Definition file not found: %s", DEF_FILE)
        return False
    if SIF_PATH.exists() and not rebuild:
        logging.info("Using cached SIF: %s", SIF_PATH)
        return True

    logging.info("Building %s from %s …", SIF_PATH, DEF_FILE)
    cmd = [SING_BIN, "build", str(SIF_PATH), str(DEF_FILE)]
    try:
        run(cmd)
        logging.info("Build finished ✓")
        return True
    except subprocess.CalledProcessError as e:
        logging.error("Singularity build failed (%s)", e.returncode)
        return False


def instance_running() -> bool:
    try:
        out = run([SING_BIN, "instance", "list"], capture=True).stdout
        return INSTANCE_NAME in out
    except subprocess.CalledProcessError:
        return False


def start_instance(rebuild: bool = False) -> bool:
    if instance_running():
        logging.warning("Instance already running – restarting…")
        stop_instance()

    if not build_sif(rebuild=rebuild):
        return False

    logging.info("Starting instance %s …", INSTANCE_NAME)
    cmd = [
        SING_BIN, "instance", "start",
        "--cleanenv",
        "--net",
        "--network-args", f"portmap={API_PORT_HOST}:tcp:{API_PORT_INSIDE}",
        str(SIF_PATH),
        INSTANCE_NAME,
    ]
    try:
        run(cmd)
        time.sleep(3)
        if instance_running():
            logging.info("Instance running. Access API at http://localhost:%d", API_PORT_HOST)
            return True
        logging.error("Instance failed to appear in list.")
        return False
    except subprocess.CalledProcessError as e:
        logging.error("Failed to start instance: %s", e)
        return False


def stop_instance() -> bool:
    if not instance_running():
        logging.info("Instance not running.")
        return True
    logging.info("Stopping instance %s …", INSTANCE_NAME)
    try:
        run([SING_BIN, "instance", "stop", INSTANCE_NAME])
        return True
    except subprocess.CalledProcessError as e:
        logging.error("Failed to stop instance: %s", e)
        return False


def show_status():
    logging.info("Instance: %s", "running" if instance_running() else "stopped")
    logging.info("API port (host): %d", API_PORT_HOST)


def show_logs(lines: int = 50):
    if not instance_running():
        logging.warning("Instance not running.")
        return
    log_dir = Path.home() / ".apptainer" / "instances" / "logs" / os.getenv("USER", "")
    log_file = log_dir / f"{INSTANCE_NAME}.log"
    if not log_file.exists():
        logging.warning("Log file not found: %s", log_file)
        return
    print("\n--- logs ---")
    print(run(["tail", "-n", str(lines), str(log_file)], capture=True).stdout)
    print("------------")

# ---------------------------------------------------------------------------
# Interactive REPL
# ---------------------------------------------------------------------------

def repl():
    print("Singularity Sandbox Manager (type 'help')")
    while True:
        try:
            line = input("cmd> ").strip()
        except EOFError:
            break
        if not line:
            continue
        cmd, *args = shlex.split(line)
        if cmd in {"exit", "quit"}:
            break
        elif cmd == "help":
            print("build | start [--rebuild] | stop | status | logs [N] | exit")
        elif cmd == "build":
            rebuild = "--rebuild" in args
            build_sif(rebuild=rebuild)
        elif cmd == "start":
            rebuild = "--rebuild" in args
            start_instance(rebuild=rebuild)
        elif cmd == "stop":
            stop_instance()
        elif cmd == "status":
            show_status()
        elif cmd == "logs":
            n = int(args[0]) if args else 50
            show_logs(n)
        else:
            print("Unknown command.")
    stop_instance()

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if len(sys.argv) == 1:
        repl()
        sys.exit(0)

    p = argparse.ArgumentParser("Singularity Sandbox Manager")
    sp = p.add_subparsers(dest="cmd", required=True)

    sp.add_parser("build").add_argument("--rebuild", action="store_true")
    sp.add_parser("start").add_argument("--rebuild", action="store_true")
    sp.add_parser("stop")
    sp.add_parser("status")
    lp = sp.add_parser("logs")
    lp.add_argument("n", nargs="?", type=int, default=50)

    a = p.parse_args()
    ok = True
    if a.cmd == "build":
        ok = build_sif(rebuild=a.rebuild)
    elif a.cmd == "start":
        ok = start_instance(rebuild=a.rebuild)
    elif a.cmd == "stop":
        ok = stop_instance()
    elif a.cmd == "status":
        show_status()
    elif a.cmd == "logs":
        show_logs(a.n)
    sys.exit(0 if ok else 1)
