#!/usr/bin/env python3
"""Singularity Sandbox Manager (Docker‑free)
==========================================
Pure‑Singularity version that **does not require Docker at all**.  It pulls a
pre-built `sandbox.sif` from a specified URL when necessary.

Commands
-------------------------
    build   – download/update `sandbox.sif` from a predefined URL
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
SIF_PATH = SCRIPT_DIR / "sandbox.sif"          # output image (will be downloaded here)
# DEF_FILE is no longer needed as we are pulling a pre-built image.
CONTAINER_URL = "https://github.com/OpenTechBio/Olaf/releases/download/v0.0.1/benchmarking_sandbox.sif"

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


def pull_sif_if_needed(force_pull: bool = False) -> bool:
    """Pull sandbox.sif from the predefined URL if it doesn't exist or if force_pull is true."""
    if SIF_PATH.exists() and not force_pull:
        logging.info("Using existing SIF: %s", SIF_PATH)
        return True

    if SIF_PATH.exists() and force_pull:
        logging.info("Force pull enabled, removing existing SIF: %s", SIF_PATH)
        try:
            SIF_PATH.unlink()
        except OSError as e:
            logging.error("Failed to remove existing SIF %s: %s", SIF_PATH, e)
            return False

    logging.info("Pulling %s from %s …", SIF_PATH, CONTAINER_URL)
    # Command: singularity pull [--force] <local_sif_name> <remote_uri>
    # Using --force to handle cases where the file might exist despite prior checks or incomplete downloads.
    cmd = [SING_BIN, "pull", "--force", str(SIF_PATH), CONTAINER_URL]
    try:
        run(cmd)
        if not SIF_PATH.exists() or SIF_PATH.stat().st_size == 0:
            logging.error("Singularity pull command executed, but SIF file is missing or empty.")
            return False
        logging.info("Pull finished ✓. SIF is at %s", SIF_PATH)
        return True
    except subprocess.CalledProcessError as e:
        logging.error("Singularity pull failed (return code %s)", e.returncode)
        # subprocess.run with text=True should populate stdout/stderr
        if hasattr(e, 'stderr') and e.stderr:
            logging.error("Stderr:\n%s", e.stderr.strip())
        if hasattr(e, 'stdout') and e.stdout: # Sometimes singularity puts error info in stdout
             logging.error("Stdout:\n%s", e.stdout.strip())
        return False
    except Exception as e:
        logging.error("An unexpected error occurred during Singularity pull: %s", e)
        return False


def instance_running() -> bool:
    try:
        # Use check=False as a non-zero exit code (no instances running) is not an error here.
        result = run([SING_BIN, "instance", "list"], capture=True, check=False)
        return INSTANCE_NAME in result.stdout
    except subprocess.CalledProcessError: # Should not happen with check=False
        return False
    except FileNotFoundError: # If SING_BIN itself is somehow removed mid-script
        logging.error("Singularity/Apptainer executable not found.")
        return False


def start_instance(rebuild: bool = False) -> bool: # rebuild here means force_pull for the SIF
    if instance_running():
        logging.warning("Instance already running – restarting…")
        if not stop_instance(): # Attempt to stop, if it fails, don't proceed.
             logging.error("Failed to stop existing instance. Cannot start new one.")
             return False

    # The 'rebuild' flag for start is interpreted as 'force_pull' for the SIF image
    if not pull_sif_if_needed(force_pull=rebuild):
        logging.error("Failed to ensure SIF image is available. Cannot start instance.")
        return False

    logging.info("Starting instance %s from %s …", INSTANCE_NAME, SIF_PATH)
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
        logging.info("Instance start command executed. Waiting a moment to check status...")
        time.sleep(3) # Give the instance a moment to register
        if instance_running():
            logging.info("Instance '%s' is running. Access API at http://localhost:%d", INSTANCE_NAME, API_PORT_HOST)
            return True
        else:
            logging.error("Instance '%s' failed to appear in list after start command.", INSTANCE_NAME)
            # Attempt to get logs if possible (might not exist if start failed very early)
            log_dir_base = Path.home() / ".apptainer" / "instances" / "logs"
            # User might not be available via os.getenv reliably in all contexts, try to find log
            # This part is heuristic for finding the log file.
            user_name = os.getenv("USER", "unknown_user")
            specific_log_dir_apptainer = log_dir_base / os.getenv("USER", "") # Apptainer specific
            specific_log_dir_singularity = Path.home() / ".singularity" / "instances" / "logs" / os.getenv("USER", "") # Older Singularity

            log_file_apptainer = specific_log_dir_apptainer / f"{INSTANCE_NAME}.log"
            log_file_singularity = specific_log_dir_singularity / f"{INSTANCE_NAME}.log"

            actual_log_file = None
            if log_file_apptainer.exists():
                actual_log_file = log_file_apptainer
            elif log_file_singularity.exists():
                 actual_log_file = log_file_singularity

            if actual_log_file:
                logging.error("Check instance logs for details: %s", actual_log_file)
                try:
                    with open(actual_log_file, "r") as lf:
                        log_tail = "".join(list(lf)[-20:]) # last 20 lines
                    logging.error("Last few log lines:\n%s", log_tail)
                except Exception as log_e:
                    logging.error("Could not read log file: %s", log_e)
            else:
                logging.error("Instance log file not found in typical locations.")
            return False
    except subprocess.CalledProcessError as e:
        logging.error("Failed to start instance (command error): %s", e)
        if hasattr(e, 'stderr') and e.stderr: logging.error("Stderr:\n%s", e.stderr.strip())
        if hasattr(e, 'stdout') and e.stdout: logging.error("Stdout:\n%s", e.stdout.strip())
        return False
    except Exception as e:
        logging.error("An unexpected error occurred trying to start the instance: %s", e)
        return False


def stop_instance() -> bool:
    if not instance_running():
        logging.info("Instance '%s' not running.", INSTANCE_NAME)
        return True
    logging.info("Stopping instance %s …", INSTANCE_NAME)
    try:
        run([SING_BIN, "instance", "stop", INSTANCE_NAME])
        # Verify it's stopped
        if not instance_running():
            logging.info("Instance '%s' stopped successfully.", INSTANCE_NAME)
            return True
        else:
            logging.warning("Instance stop command executed, but instance still appears in list. Check manually.")
            return False # Or True depending on desired strictness
    except subprocess.CalledProcessError as e:
        logging.error("Failed to stop instance: %s", e)
        if hasattr(e, 'stderr') and e.stderr: logging.error("Stderr:\n%s", e.stderr.strip())
        if hasattr(e, 'stdout') and e.stdout: logging.error("Stdout:\n%s", e.stdout.strip())
        return False
    except Exception as e:
        logging.error("An unexpected error occurred trying to stop the instance: %s", e)
        return False


def show_status():
    is_running = instance_running()
    logging.info("Instance '%s': %s", INSTANCE_NAME, "RUNNING" if is_running else "STOPPED")
    if is_running:
        logging.info("API access (host): http://localhost:%d", API_PORT_HOST)
    else:
        logging.info("API port (host - if running): %d", API_PORT_HOST)


def show_logs(lines: int = 50):
    # Determine log directory based on Apptainer/Singularity conventions
    # Singularity: $HOME/.singularity/instances/logs/<USER>/<INSTANCE_NAME>.log
    # Apptainer:   $HOME/.apptainer/instances/logs/<USER>/<INSTANCE_NAME>.log
    # Prefer Apptainer path if it exists, fallback to Singularity
    user_name = os.getenv("USER")
    if not user_name:
        logging.error("USER environment variable not set, cannot reliably determine log path.")
        # Fallback for some systems where USER might not be set in certain execution contexts
        try:
            import getpass
            user_name = getpass.getuser()
        except Exception:
            logging.error("Could not determine username to find logs.")
            return

    log_dir_apptainer = Path.home() / ".apptainer" / "instances" / "logs" / user_name
    log_dir_singularity = Path.home() / ".singularity" / "instances" / "logs" / user_name

    log_file_apptainer = log_dir_apptainer / f"{INSTANCE_NAME}.log"
    log_file_singularity = log_dir_singularity / f"{INSTANCE_NAME}.log"

    actual_log_file = None
    if log_file_apptainer.exists():
        actual_log_file = log_file_apptainer
    elif log_file_singularity.exists():
        actual_log_file = log_file_singularity

    if not instance_running() and not actual_log_file : # If not running and no log file, nothing to show
        logging.warning("Instance '%s' not running and no log file found.", INSTANCE_NAME)
        return
    elif not actual_log_file: # Running but somehow no log file yet (or path issue)
        logging.warning("Instance '%s' is running, but its log file was not found at expected locations:\n- %s\n- %s", INSTANCE_NAME, log_file_apptainer, log_file_singularity)
        return
    elif not instance_running() and actual_log_file:
        logging.info("Instance '%s' is not running. Displaying last logs from %s:", INSTANCE_NAME, actual_log_file)


    print(f"\n--- Logs for {INSTANCE_NAME} (last {lines} lines from {actual_log_file}) ---")
    try:
        log_content = run(["tail", "-n", str(lines), str(actual_log_file)], capture=True, check=True).stdout
        print(log_content.strip())
    except subprocess.CalledProcessError as e:
        logging.error("Could not read logs using tail: %s", e)
        if hasattr(e, 'stderr') and e.stderr: logging.error("Stderr:\n%s", e.stderr.strip())
    except FileNotFoundError: # if tail is not found
        logging.error("'tail' command not found. Cannot display logs.")
    print("-------------------------------------")

# ---------------------------------------------------------------------------
# Interactive REPL
# ---------------------------------------------------------------------------

def repl():
    print("Singularity Sandbox Manager (pulls pre-built SIF). Type 'help'.")
    while True:
        try:
            line = input(f"{INSTANCE_NAME}> ").strip()
        except EOFError:
            print("\nExiting.")
            break
        if not line:
            continue
        
        try:
            cmd_parts = shlex.split(line)
            cmd = cmd_parts[0]
            args = cmd_parts[1:]
        except ValueError as e: # Handle issue with shlex.split if quotes are mismatched
            print(f"Error parsing command: {e}")
            continue

        if cmd in {"exit", "quit", "q"}:
            print("Exiting.")
            break
        elif cmd == "help":
            print("\nAvailable commands:")
            print("  build [--rebuild]     - Ensure SIF image is downloaded (use --rebuild to force re-download).")
            print("  start [--rebuild]     - Start the instance (forces SIF re-download if --rebuild is used).")
            print("  stop                  - Stop & remove the instance.")
            print("  status                - Show instance + port status.")
            print("  logs [N]              - Tail the instance log (default 50 lines).")
            print("  exit | quit | q       - Exit the manager.\n")
        elif cmd == "build":
            rebuild_flag = "--rebuild" in args
            pull_sif_if_needed(force_pull=rebuild_flag)
        elif cmd == "start":
            rebuild_flag = "--rebuild" in args
            start_instance(rebuild=rebuild_flag)
        elif cmd == "stop":
            stop_instance()
        elif cmd == "status":
            show_status()
        elif cmd == "logs":
            n_lines = 50
            if args:
                try:
                    n_lines = int(args[0])
                except ValueError:
                    print("Invalid number of lines. Using default 50.")
            show_logs(n_lines)
        else:
            print(f"Unknown command: {cmd}. Type 'help' for available commands.")
    
    # Attempt to gracefully stop the instance if it's running when REPL exits
    if instance_running():
        logging.info("REPL exited, stopping instance '%s' if running...", INSTANCE_NAME)
        stop_instance()

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Configure logging
    # Add a timestamp to the logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")


    if len(sys.argv) == 1: # No arguments, start REPL
        repl()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="Singularity Sandbox Manager. Pulls a pre-built SIF and manages its instance.",
        formatter_class=argparse.RawTextHelpFormatter # To preserve help text formatting
    )
    subparsers = parser.add_subparsers(dest="cmd", title="commands", required=True)

    # Build command (now means: ensure SIF is downloaded)
    build_parser = subparsers.add_parser("build", help="Ensure SIF image is downloaded from the predefined URL.")
    build_parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Force re-download of the SIF image even if it exists."
    )

    # Start command
    start_parser = subparsers.add_parser("start", help="Start the Singularity instance.")
    start_parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Force re-download of the SIF image before starting."
    )

    # Stop command
    stop_parser = subparsers.add_parser("stop", help="Stop and remove the Singularity instance.")

    # Status command
    status_parser = subparsers.add_parser("status", help="Show current status of the instance and API port.")

    # Logs command
    logs_parser = subparsers.add_parser("logs", help="Show logs from the running instance.")
    logs_parser.add_argument(
        "n",
        nargs="?",
        type=int,
        default=50,
        help="Number of log lines to display (default: 50)."
    )

    args = parser.parse_args()
    operation_successful = True # Assume success unless a command returns False

    if args.cmd == "build":
        operation_successful = pull_sif_if_needed(force_pull=args.rebuild)
    elif args.cmd == "start":
        operation_successful = start_instance(rebuild=args.rebuild)
    elif args.cmd == "stop":
        operation_successful = stop_instance()
    elif args.cmd == "status":
        show_status() # Status typically doesn't fail in a way that sets operation_successful
    elif args.cmd == "logs":
        show_logs(args.n) # Logs display also doesn't typically set operation_successful

    sys.exit(0 if operation_successful else 1)