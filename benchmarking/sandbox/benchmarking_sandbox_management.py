import docker
import argparse
import os
import sys
import time
import subprocess # Still needed for docker cp
import shlex
import json
import tarfile
import io
import tempfile # Needed for temporary connection file

# --- Dependency Imports ---
try:
    import jupyter_client
    from jupyter_client.blocking.client import BlockingKernelClient # Explicit import
    from queue import Empty # For checking kernel message queue
    from jupyter_client.kernelspec import KernelSpecManager # To check for kernelspec if needed
except ImportError:
    print("Error: jupyter_client library not found.", file=sys.stderr)
    print("Please install it in your host environment: pip install jupyter_client", file=sys.stderr)
    sys.exit(1)

# Optional: Use rich for better formatting
try:
    from rich.console import Console
    from rich.prompt import Prompt
    HAS_RICH = True
except ImportError:
    HAS_RICH = False
    class Console:
        def print(self, *args, **kwargs):
             kwargs.pop('style', None); kwargs.pop('justify', None)
             if 'file' in kwargs: print(*args, **kwargs)
             else: print(*args)
    class Prompt:
        @staticmethod
        def ask(prompt, choices=None, default=None):
            p_text = f"{prompt} "
            if choices: p_text += f"({'/'.join(choices)}) "
            if default: p_text += f"[{default}] "
            return input(p_text).strip()

# --- Configuration ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCKERFILE_PATH = os.path.join(SCRIPT_DIR, 'Dockerfile')
IMAGE_TAG = "benchmarking-sandbox:latest"
CONTAINER_NAME = "benchmarking_sandbox_instance"
KERNEL_CONNECTION_FILE_PATH = "/home/sandboxuser/kernel-connection.json" # Path inside container
FIXED_BASE_PORT = 4000          # keep in sync with Dockerfile

# Instantiate console globally
console = Console() if HAS_RICH else None

def _print_message(message, style=None, is_error=False):
     """Helper to print using rich console or standard print, handling stderr."""
     if is_error: print(message, file=sys.stderr)
     elif console: console.print(message, style=style)
     else: print(message)

class SandboxManager:
    """
    Manages the lifecycle and kernel communication of the benchmarking sandbox Docker container.
    """
    def __init__(self):
        self.client = None
        self.container = None
        self.kernel_manager = None
        self.kernel_client = None
        self.connection_info = None
        self._temp_conn_file_path = None
        try:
            docker_host = os.environ.get("DOCKER_HOST")
            if docker_host: _print_message(f"Using DOCKER_HOST: {docker_host}")
            self.client = docker.from_env()
            self.client.ping()
            _print_message("Docker client initialized successfully.")
        except Exception as e:
            _print_message(f"Error initializing Docker client: {e}", is_error=True)
            _print_message("Ensure Docker Desktop/Engine is running and DOCKER_HOST is set if needed.", is_error=True)
            sys.exit(1)

    def _get_container_logs(self, tail=50):
        """Retrieves logs from the container."""
        if not self.container:
            return "(Container object not available)"
        try:
            self.container.reload() # Get latest state
            logs = self.container.logs(tail=tail).decode('utf-8', errors='ignore')
            return logs
        except docker.errors.NotFound:
             return "(Container not found for logs)"
        except Exception as log_e:
             return f"(Could not retrieve logs: {log_e})"

    def _find_container(self):
        """Finds the container by name, returns container object or None."""
        try:
            # Don't store self.container here, let other methods do that
            container = self.client.containers.get(CONTAINER_NAME)
            return container
        except docker.errors.NotFound:
            return None
        except Exception as e:
            _print_message(f"Error finding container '{CONTAINER_NAME}': {e}", is_error=True)
            return None

    def build_image(self):
        """Builds the Docker image from the Dockerfile."""
        # (Build logic remains largely the same)
        _print_message(f"Building Docker image '[cyan]{IMAGE_TAG}[/cyan]' from [blue]{DOCKERFILE_PATH}[/blue]...", style="bold")
        if not os.path.exists(DOCKERFILE_PATH):
            _print_message(f"Error: Dockerfile not found at {DOCKERFILE_PATH}", style="bold red", is_error=True); return False
        try:
            build_context = os.path.dirname(DOCKERFILE_PATH)
            _print_message(f"Using build context: [blue]{build_context}[/blue]")
            stream = self.client.api.build(path=build_context, dockerfile=os.path.basename(DOCKERFILE_PATH), tag=IMAGE_TAG, rm=True, decode=True)
            last_status = None
            for chunk in stream:
                if 'stream' in chunk:
                    line = chunk['stream'].strip();
                    if line: _print_message(line)
                elif 'errorDetail' in chunk:
                    _print_message(f"Build Error: {chunk['errorDetail']['message']}", style="bold red", is_error=True); return False
                elif 'status' in chunk:
                    status = chunk['status']
                    if status != last_status: _print_message(f"Status: {status}"); last_status = status
            _print_message(f"Image '[cyan]{IMAGE_TAG}[/cyan]' built successfully.", style="green")
            return True
        except docker.errors.BuildError as e:
            _print_message(f"Docker build failed: {e}", style="bold red", is_error=True)
            for line in e.build_log:
                 if 'stream' in line: _print_message(line['stream'].strip(), style="dim", is_error=True)
            return False
        except Exception as e:
            _print_message(f"An unexpected error occurred during build: {e}", style="bold red", is_error=True); return False

    def start_container(self, rebuild=False):
            """Starts the Docker container which runs the kernel and retrieves connection info."""
            if self.kernel_client and self.kernel_client.is_alive():
                _print_message("Kernel connection already active.", style="yellow"); return True

            if rebuild:
                _print_message("Rebuild requested.")
                existing_container = self._find_container()
                if existing_container:
                    _print_message("Stopping existing container before rebuild...")
                    self.stop_container(remove=True, container_obj=existing_container)
                if not self.build_image():
                    _print_message("Build failed, cannot start container.", style="bold red", is_error=True); return False

            existing_container = self._find_container()
            if existing_container:
                _print_message(f"Found existing container '{CONTAINER_NAME}'. Stopping and removing it first...")
                if not self.stop_container(remove=True, container_obj=existing_container):
                    _print_message(f"Failed to stop/remove existing container '{CONTAINER_NAME}'. Aborting.", style="bold red", is_error=True); return False

            _print_message(f"Starting container '[cyan]{CONTAINER_NAME}[/cyan]' to run kernel...")
            try:
                try: self.client.images.get(IMAGE_TAG)
                except docker.errors.ImageNotFound:
                    _print_message(f"Image '[cyan]{IMAGE_TAG}[/cyan]' not found. Building...", style="yellow")
                    if not self.build_image(): _print_message("Build failed.", style="bold red", is_error=True); return False
                # Add in port mapping
                port_map = {f"{p}/tcp": FIXED_BASE_PORT + i
                for i, p in enumerate(range(FIXED_BASE_PORT, FIXED_BASE_PORT + 5))}
                # Store container object in the instance variable
                self.container = self.client.containers.run(IMAGE_TAG, name=CONTAINER_NAME, detach=True, auto_remove=False, ports=port_map)
                # --------> REMOVED THE ERRONEOUS LINE FROM HERE <----------
                # self.connection_info = json.loads(connection_file_content) # <<-- REMOVED

                _print_message(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' started ([yellow]{self.container.short_id}[/yellow]).")
                _print_message("Waiting for kernel connection file...")

                max_wait = 30; start_time = time.time()
                connection_file_content = None # Initialize to None before the loop
                while time.time() - start_time < max_wait:
                    try:
                        # Check container status before getting archive
                        self.container.reload()
                        if self.container.status != 'running':
                            _print_message(f"Container exited unexpectedly (status: {self.container.status}).", style="bold red", is_error=True)
                            logs = self._get_container_logs()
                            _print_message("--- Container Logs ---", is_error=True)
                            _print_message(logs, is_error=True)
                            _print_message("----------------------", is_error=True)
                            # Attempt cleanup even on unexpected exit before returning False
                            if self.container: self.stop_container(remove=True, container_obj=self.container)
                            return False

                        bits, stat = self.container.get_archive(KERNEL_CONNECTION_FILE_PATH)
                        file_obj = io.BytesIO();
                        for chunk in bits: file_obj.write(chunk)
                        file_obj.seek(0)
                        with tarfile.open(fileobj=file_obj, mode='r') as tar:
                            member = tar.getmembers()[0]
                            extracted_file = tar.extractfile(member)
                            if extracted_file:
                                connection_file_content = extracted_file.read().decode('utf-8') # Content assigned here
                                _print_message("[green]Kernel connection file retrieved.[/green]")
                                break # Exit loop on success
                    except docker.errors.NotFound:
                        # File not ready yet, container still running, wait a bit
                        time.sleep(0.5)
                    except Exception as e:
                        _print_message(f"Error retrieving connection file: {e}", style="red", is_error=True)
                        # Attempt cleanup before returning False
                        if self.container: self.stop_container(remove=True, container_obj=self.container)
                        return False
                # This 'else' block executes if the while loop finishes WITHOUT hitting 'break'
                else:
                    _print_message(f"Kernel connection file not found after {max_wait} seconds.", style="bold red", is_error=True)
                    logs = self._get_container_logs()
                    _print_message("--- Container Logs ---", is_error=True)
                    _print_message(logs, is_error=True)
                    _print_message("----------------------", is_error=True)
                    # Attempt cleanup before returning False
                    if self.container: self.stop_container(remove=True, container_obj=self.container)
                    return False

                # --------> CORRECT PLACE TO PARSE THE JSON <----------
                # Ensure content was actually retrieved before trying to parse
                if connection_file_content:
                    self.connection_info = json.loads(connection_file_content)
                    _print_message(f"Kernel ports: "
                    f"{[self.connection_info[p] for p in ('shell_port','iopub_port','stdin_port','hb_port','control_port')]}")
                    # self.connection_info["ip"] = "127.0.0.1" # Override IP for localhost access
                    # Override ports based on the mapping used when running the container
                    # for off, key in enumerate(
                    #         ["shell_port", "iopub_port", "stdin_port", "hb_port", "control_port"]):
                    #     self.connection_info[key] = FIXED_BASE_PORT + off
                    return True # Success!
                else:
                    # This case should ideally not be reached if the loop logic is correct,
                    # but handle defensively.
                    _print_message("Failed to retrieve connection file content.", style="bold red", is_error=True)
                    if self.container: self.stop_container(remove=True, container_obj=self.container)
                    return False


            except Exception as e:
                _print_message(f"Error starting container/kernel: {e}", style="bold red", is_error=True)
                # Ensure cleanup if self.container was assigned during the run attempt
                current_container = self._find_container() # Re-check by name in case self.container wasn't set or is stale
                if current_container:
                    self.stop_container(remove=True, container_obj=current_container)
                else: # If container object isn't available, clear internal state
                    self.container = None
                return False

    def connect_kernel(self):
        """Connects to the running kernel using the retrieved connection info."""
        if not self.connection_info:
             _print_message("Cannot connect: No kernel connection info available.", style="red", is_error=True); return False
        if self.kernel_client and self.kernel_client.is_alive():
             _print_message("Already connected to kernel.", style="yellow"); return True

        temp_conn_file = None
        try:
            _print_message("Connecting to kernel...")
            temp_conn_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json', encoding='utf-8')
            json.dump(self.connection_info, temp_conn_file)
            temp_conn_file.close()
            self._temp_conn_file_path = temp_conn_file.name
            _print_message(f"Wrote connection info to temporary file: {self._temp_conn_file_path}", style="dim")

            self.kernel_client = BlockingKernelClient(connection_file=self._temp_conn_file_path)
            self.kernel_client.start_channels()
            _print_message("Kernel client channels started.")

            # Test connection - Increased timeout, added log retrieval on failure
            connection_timeout = 20
            _print_message(f"Waiting for kernel ready (timeout={connection_timeout}s)...")
            self.kernel_client.wait_for_ready(timeout=connection_timeout)
            _print_message("[green]Successfully connected to kernel.[/green]")
            return True
        # *** MODIFICATION: Catch specific RuntimeError for kernel death ***
        except RuntimeError as e:
             if "Kernel died" in str(e):
                  _print_message(f"Failed to connect to kernel: {e}", style="bold red", is_error=True)
                  logs = self._get_container_logs()
                  _print_message("--- Container Logs ---", is_error=True)
                  _print_message(logs, is_error=True)
                  _print_message("----------------------", is_error=True)
             else: # Re-raise other RuntimeErrors
                  _print_message(f"Runtime error connecting to kernel: {e}", style="bold red", is_error=True)
             self.kernel_client = None
             self._cleanup_temp_conn_file()
             return False
        except Exception as e:
            _print_message(f"Failed to connect to kernel: {e}", style="bold red", is_error=True)
            # Attempt to get logs if container object exists
            if self.container:
                 logs = self._get_container_logs()
                 _print_message("--- Container Logs ---", is_error=True)
                 _print_message(logs, is_error=True)
                 _print_message("----------------------", is_error=True)
            self.kernel_client = None
            self._cleanup_temp_conn_file()
            return False
        # *********************************************************

    def _cleanup_temp_conn_file(self):
         """Cleans up the temporary connection file if it exists."""
         # (Cleanup logic remains the same)
         if self._temp_conn_file_path and os.path.exists(self._temp_conn_file_path):
              try:
                   os.remove(self._temp_conn_file_path)
                   _print_message(f"Cleaned up temporary connection file: {self._temp_conn_file_path}", style="dim")
                   self._temp_conn_file_path = None
              except OSError as e:
                   _print_message(f"Warning: Could not remove temp file {self._temp_conn_file_path}: {e}", style="yellow", is_error=True)

    def disconnect_kernel(self):
        """Disconnects from the kernel and cleans up temp file."""
        # (Disconnect logic remains the same)
        if self.kernel_client:
            _print_message("Disconnecting from kernel...")
            try:
                if self.kernel_client.is_alive(): self.kernel_client.stop_channels()
                _print_message("Kernel client channels stopped.")
            except Exception as e: _print_message(f"Error stopping kernel channels: {e}", style="yellow", is_error=True)
            self.kernel_client = None; self.connection_info = None
        else: _print_message("Not connected to any kernel.")
        self._cleanup_temp_conn_file()


    def stop_container(self, remove=False, container_obj=None):
        return
        """Stops the container and disconnects the kernel."""
        # (Stop logic remains largely the same, ensures disconnect first)
        self.disconnect_kernel()
        container_to_stop = container_obj or self._find_container()
        if not container_to_stop:
            _print_message(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' not found.")
            self.container = None # Ensure cleared
            return True
        _print_message(f"Stopping container '[cyan]{CONTAINER_NAME}[/cyan]' ([yellow]{container_to_stop.short_id}[/yellow])...")
        stopped = False; removed_flag = False
        try:
            current_status = container_to_stop.status
            if current_status == 'running':
                container_to_stop.stop(timeout=5); _print_message("Container stop signal sent.")
                container_to_stop.wait(timeout=10); _print_message("Container stopped.")
            else: _print_message(f"Container was not running (status: {current_status}).")
            stopped = True
            if remove:
                try:
                     container_to_stop.remove(force=True); _print_message("Container removed.")
                     removed_flag = True
                except docker.errors.APIError as e:
                     if e.response.status_code != 404: _print_message(f"API error removing: {e}", style="yellow", is_error=True)
                     else: _print_message("Container was already removed.")
                     removed_flag = True
            else: removed_flag = True
            self.container = None
            return stopped and removed_flag
        except Exception as e:
            _print_message(f"Error stopping/removing container: {e}", style="bold red", is_error=True)
            if remove and container_to_stop:
                 try: container_to_stop.remove(force=True); _print_message("Force removed container.")
                 except: pass
            self.container = None
            return False

    def get_status(self):
        """Gets the status of the container and kernel connection."""
        # (Status logic remains largely the same)
        container_status = "not found"; kernel_status = "disconnected"
        container = self._find_container() # Use local var, don't rely on self.container always being set
        if container:
            try: container.reload(); container_status = container.status
            except Exception: container_status = "unknown (error)"
        if self.kernel_client and self.kernel_client.is_alive(): kernel_status = "connected"
        elif self.connection_info and container_status == 'running':
             try: # Check if kernel process is running inside
                  exit_code, output = container.exec_run("pgrep -f ipykernel_launcher")
                  if exit_code == 0: kernel_status = "pending connection"
                  else: kernel_status = "kernel process not found"
             except Exception: kernel_status = "unknown (exec check failed)"
        return f"Container: {container_status}, Kernel: {kernel_status}"


    def run_code(self, code_string, timeout=60):
        """Runs a Python code string using the connected Jupyter kernel."""
        # (Run code logic remains the same)
        if not self.kernel_client or not self.kernel_client.is_alive():
            _print_message("Error: Not connected to a running kernel.", style="bold red", is_error=True); return None
        _print_message(f"Executing code via kernel '[cyan]{self.kernel_client.kernel_id}[/cyan]'...")
        msg_id = self.kernel_client.execute(code_string, store_history=True)
        stdout_list = []; stderr_list = []; finished = False; start_time = time.time()
        while time.time() - start_time < timeout:
            try:
                msg = self.kernel_client.get_iopub_msg(timeout=1)
                msg_type = msg['header']['msg_type']; content = msg['content']
                if msg_type == 'status':
                    if content['execution_state'] == 'idle': finished = True; break
                elif msg_type == 'stream':
                    if content['name'] == 'stdout': stdout_list.append(content['text'])
                    elif content['name'] == 'stderr': stderr_list.append(content['text'])
                elif msg_type == 'execute_result':
                    if 'data' in content and 'text/plain' in content['data']: stdout_list.append(content['data']['text/plain'] + '\n')
                elif msg_type == 'display_data':
                     if 'data' in content and 'text/plain' in content['data']: stdout_list.append(content['data']['text/plain'] + '\n')
                elif msg_type == 'error':
                    tb = content.get('traceback', [])
                    stderr_list.append(f"Error: {content.get('ename', 'Err')}: {content.get('evalue', '')}\n")
                    stderr_list.extend(line + '\n' for line in tb); finished = True; break
            except Empty: pass
            except Exception as e: _print_message(f"Error reading kernel message: {e}", style="red", is_error=True); break
            if finished: break
        final_status = "unknown"
        if finished:
             try:
                  reply = self.kernel_client.get_shell_msg(timeout=5)
                  if reply['header']['msg_type'] == 'execute_reply': final_status = reply['content']['status']
                  else: _print_message(f"Unexpected final reply type: {reply['header']['msg_type']}", style="yellow")
             except Empty: _print_message("Timeout or no reply on shell channel after execution.", style="yellow", is_error=True)
             except Exception as e: _print_message(f"Error reading final kernel reply: {e}", style="red", is_error=True)
        elif time.time() - start_time >= timeout: _print_message(f"Code execution timed out after {timeout} seconds.", style="bold red", is_error=True); return None
        final_stdout = "".join(stdout_list).strip(); final_stderr = "".join(stderr_list).strip()
        _print_message("[bold]--- Execution Output (stdout) ---[/bold]")
        if final_stdout: _print_message(final_stdout)
        else: _print_message("[No standard output captured]")
        _print_message("[bold]---------------------------------[/bold]")
        if final_stderr:
             _print_message("[bold red]--- Execution Errors (stderr) ---[/bold red]", is_error=True)
             _print_message(final_stderr, is_error=True)
             _print_message("[bold red]---------------------------------[/bold red]", is_error=True)
        return final_stdout


# --- Interactive Mode Functions ---

def print_interactive_help():
    """Prints help message for interactive mode."""
    # (Help text remains the same)
    _print_message("\n[bold cyan]Available Commands:[/bold cyan]")
    _print_message("  [green]build[/green]               Build the Docker image.")
    _print_message("  [green]start[/green] [--rebuild]   Start container & kernel (optionally rebuild).")
    _print_message("  [green]connect[/green]             Connect to kernel in running container.")
    _print_message("  [green]stop[/green]                Disconnect kernel, stop & remove container.")
    _print_message("  [green]status[/green]              Check container and kernel status.")
    _print_message("  [green]run[/green] <code>          Run Python code via connected kernel.")
    _print_message("                      (Example: run \"print('hello')\")")
    _print_message("  [green]help[/green]                Show this help message.")
    _print_message("  [green]exit[/green]                Exit (prompts to stop container if running).")
    _print_message("\nExample: [yellow]start --rebuild[/yellow]")
    _print_message("Example: [yellow]run \"import numpy as np; a = np.array([1,2]); print(a*2)\"[/yellow]")

def interactive_loop(manager):
    """Runs the interactive command loop."""
    # (Interactive loop logic remains the same)
    _print_message("[bold blue]Welcome to the Stateful Benchmarking Sandbox Manager![/bold blue]")
    print_interactive_help()
    while True:
        try:
            if HAS_RICH: raw_command = Prompt.ask("\nEnter command (\'help\' or \'exit\')")
            else: raw_command = input("\nEnter command ('help' or 'exit'): ").strip()
            if not raw_command: continue
            try: command_parts = shlex.split(raw_command)
            except ValueError as e: _print_message(f"Error parsing command: {e}", style="red", is_error=True); continue
            if not command_parts: continue
            command = command_parts[0].lower(); args = command_parts[1:]
            if command == "exit":
                container_running = False; container_obj = manager._find_container()
                if container_obj: 
                    try: 
                        container_obj.reload()
                        container_running = container_obj.status == 'running' 
                    except Exception:
                        pass
                if container_running:
                     should_stop = Prompt.ask(f"Container '{CONTAINER_NAME}' is running. Stop it?", choices=["y", "n"], default="y").lower() == 'y'
                     if should_stop: manager.stop_container(remove=True, container_obj=container_obj)
                break
            elif command == "help": print_interactive_help()
            elif command == "build":
                if len(args) == 0: manager.build_image()
                else: _print_message("Usage: build", style="yellow")
            elif command == "start":
                rebuild = '--rebuild' in args
                if all(a == '--rebuild' for a in args) or len(args) == 0:
                     if manager.start_container(rebuild=rebuild): manager.connect_kernel()
                else: _print_message("Usage: start [--rebuild]", style="yellow")
            elif command == "connect":
                 if len(args) == 0: manager.connect_kernel()
                 else: _print_message("Usage: connect", style="yellow")
            elif command == "stop":
                 if len(args) == 0: manager.stop_container(remove=True)
                 else: _print_message("Usage: stop", style="yellow")
            elif command == "status":
                 if len(args) == 0: _print_message(f"Status: {manager.get_status()}")
                 else: _print_message("Usage: status", style="yellow")
            elif command == "run":
                 if len(args) == 1: manager.run_code(args[0])
                 else: _print_message("Usage: run \"<python_code_string>\"", style="yellow")
            else: _print_message(f"Unknown command: '{command}'. Type 'help'.", style="red")
        except EOFError:
             _print_message("\nEOF detected. Exiting.", style="yellow")
             container_running = False; container_obj = manager._find_container()
             if container_obj: 
                try: 
                    container_obj.reload(); 
                    container_running = (container_obj.status == 'running')
                except Exception: 
                    pass
             if container_running:
                  should_stop = Prompt.ask(f"Container '{CONTAINER_NAME}' is running. Stop it?", choices=["y", "n"], default="y").lower() == 'y'
                  if should_stop: manager.stop_container(remove=True, container_obj=container_obj)
             break
        except KeyboardInterrupt: _print_message("\nInterrupted by user. Type 'exit'.", style="yellow")
        except Exception as e: _print_message(f"Unexpected error: {e}", style="bold red", is_error=True)
    _print_message("Exiting sandbox manager.", style="bold blue")


# --- Main Entry Point ---
def main():
    # (Main entry logic remains the same)
    try: manager = SandboxManager()
    except SystemExit: sys.exit(1)
    if len(sys.argv) == 1: interactive_loop(manager); sys.exit(0)
    parser = argparse.ArgumentParser(description="Manage the Stateful Benchmarking Sandbox.", formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', help='Action to perform', required=True)
    parser_build = subparsers.add_parser('build', help='Build the Docker image.')
    parser_build.set_defaults(func=lambda args, mgr: mgr.build_image())
    parser_start = subparsers.add_parser('start', help='Start container (kernel runs automatically).')
    parser_start.add_argument('--rebuild', action='store_true', help='Rebuild image first.')
    parser_start.set_defaults(func=lambda args, mgr: mgr.start_container(rebuild=args.rebuild))
    parser_connect = subparsers.add_parser('connect', help='Connect to kernel in running container.')
    parser_connect.set_defaults(func=lambda args, mgr: mgr.connect_kernel())
    parser_stop = subparsers.add_parser('stop', help='Stop container & disconnect kernel.')
    parser_stop.set_defaults(func=lambda args, mgr: mgr.stop_container(remove=True))
    parser_status = subparsers.add_parser('status', help='Check container and kernel status.')
    parser_status.set_defaults(func=lambda args, mgr: _print_message(f"Status: {mgr.get_status()}"))
    parser_run = subparsers.add_parser('run', help='Connect, run Python code via kernel, then disconnect.')
    parser_run.add_argument('code', type=str, help='The Python code string to execute.')
    def run_and_disconnect(args, mgr):
        if not mgr.connect_kernel(): return False
        result = mgr.run_code(args.code)
        mgr.disconnect_kernel()
        return result is not None
    parser_run.set_defaults(func=run_and_disconnect)
    args = parser.parse_args()
    result = args.func(args, manager)
    sys.exit(0 if result else 1)


if __name__ == "__main__":
    main()
