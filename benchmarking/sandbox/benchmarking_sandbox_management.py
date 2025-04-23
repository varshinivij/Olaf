import docker
import argparse
import os
import sys
import time
import subprocess # For running docker exec
import shlex      # For parsing interactive commands

# Optional: Use rich for better formatting
try:
    from rich.console import Console
    from rich.prompt import Prompt # For interactive prompts
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
            if choices:
                choices_str = '/'.join(choices)
                p_text += f"({choices_str}) "
            if default:
                p_text += f"[{default}] "
            return input(p_text).strip()

# --- Configuration ---
# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCKERFILE_PATH = os.path.join(SCRIPT_DIR, 'Dockerfile') # Assumes Dockerfile is in the same dir
# Use a consistent image tag and container name
IMAGE_TAG = "benchmarking-sandbox:latest"
CONTAINER_NAME = "benchmarking_sandbox_instance"
JUPYTER_PORT = 8888 # Port used internally and potentially mapped

class SandboxManager:
    """Manages the lifecycle of the benchmarking sandbox Docker container."""

    def __init__(self):
        try:
            # Check for DOCKER_HOST env var, otherwise use default
            # Note: docker.from_env() handles this automatically,
            # but explicit check can be useful for debugging.
            docker_host = os.environ.get("DOCKER_HOST")
            if docker_host:
                 print(f"Using DOCKER_HOST from environment: {docker_host}")
            self.client = docker.from_env()
            # Test connection
            self.client.ping()
            print("Docker client initialized successfully.")
        except Exception as e:
            print(f"Error initializing Docker client: {e}", file=sys.stderr)
            print("Please ensure Docker Desktop or Docker Engine is running and accessible.", file=sys.stderr)
            print("If using a non-standard Docker socket, ensure DOCKER_HOST environment variable is set.", file=sys.stderr)
            sys.exit(1)

    def _find_container(self):
        """Finds the container by name, returns container object or None."""
        try:
            return self.client.containers.get(CONTAINER_NAME)
        except docker.errors.NotFound:
            return None
        except Exception as e:
            print(f"Error finding container '{CONTAINER_NAME}': {e}", file=sys.stderr)
            return None # Treat other errors as container not found for safety

    def build_image(self):
        """Builds the Docker image from the Dockerfile."""
        console = Console() # Use Console for rich output if available
        console.print(f"Building Docker image '[cyan]{IMAGE_TAG}[/cyan]' from [blue]{DOCKERFILE_PATH}[/blue]...")
        if not os.path.exists(DOCKERFILE_PATH):
            console.print(f"[bold red]Error:[/bold red] Dockerfile not found at {DOCKERFILE_PATH}")
            return False

        try:
            build_context = os.path.dirname(DOCKERFILE_PATH)
            console.print(f"Using build context: [blue]{build_context}[/blue]")

            stream = self.client.api.build(
                path=build_context,
                dockerfile=os.path.basename(DOCKERFILE_PATH),
                tag=IMAGE_TAG,
                rm=True,
                decode=True
            )

            last_status = None
            for chunk in stream:
                if 'stream' in chunk:
                    line = chunk['stream'].strip()
                    if line: # Avoid printing empty lines
                        console.print(line)
                elif 'errorDetail' in chunk:
                    console.print(f"[bold red]Build Error:[/bold red] {chunk['errorDetail']['message']}")
                    return False
                elif 'status' in chunk:
                    status = chunk['status']
                    # Only print status if it changes to avoid clutter
                    if status != last_status:
                        console.print(f"Status: {status}")
                        last_status = status

            console.print(f"Image '[cyan]{IMAGE_TAG}[/cyan]' built successfully.")
            return True

        except docker.errors.BuildError as e:
            console.print(f"[bold red]Docker build failed:[/bold red] {e}")
            for line in e.build_log:
                 if 'stream' in line:
                     console.print(line['stream'].strip(), style="dim") # Print build log dimmed
            return False
        except docker.errors.APIError as e:
            console.print(f"[bold red]Docker API error during build:[/bold red] {e}")
            return False
        except Exception as e:
            console.print(f"[bold red]An unexpected error occurred during build:[/bold red] {e}")
            return False

    def start_container(self, rebuild=False):
        """Starts the Docker container."""
        console = Console()
        if rebuild:
            console.print("Rebuild requested.")
            if not self.build_image():
                console.print("[bold red]Build failed, cannot start container.[/bold red]")
                return None # Return None on failure

        container = self._find_container()
        if container:
            if container.status == 'running':
                console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' is already running.")
                return container
            elif container.status in ['exited', 'created']:
                console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' exists but is stopped. Starting it...")
                try:
                    container.start()
                    console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' started.")
                    console.print("Waiting a few seconds for internal services...")
                    time.sleep(5) # Wait a bit longer
                    return container
                except Exception as e:
                    console.print(f"[bold red]Error starting existing container:[/bold red] {e}")
                    return None
            else:
                console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' found in unexpected state: {container.status}. Attempting removal...")
                if not self.stop_container(): # Try stopping and removing first
                     console.print("[bold red]Failed to remove container in unexpected state. Cannot proceed.[/bold red]")
                     return None

        console.print(f"Starting a new container '[cyan]{CONTAINER_NAME}[/cyan]' from image '[cyan]{IMAGE_TAG}[/cyan]'...")
        try:
            try:
                self.client.images.get(IMAGE_TAG)
            except docker.errors.ImageNotFound:
                console.print(f"Image '[cyan]{IMAGE_TAG}[/cyan]' not found. Attempting to build...")
                if not self.build_image():
                    console.print("[bold red]Build failed, cannot start container.[/bold red]")
                    return None

            container = self.client.containers.run(
                IMAGE_TAG,
                name=CONTAINER_NAME,
                ports={f'{JUPYTER_PORT}/tcp': JUPYTER_PORT},
                detach=True,
                auto_remove=False,
            )
            console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' started with ID: [yellow]{container.short_id}[/yellow]")
            console.print(f"Waiting a few seconds for internal services (Jupyter) to potentially start...")
            time.sleep(10)
            console.print("Container should be ready.")
            return container
        except docker.errors.APIError as e:
            console.print(f"[bold red]Docker API error starting container:[/bold red] {e}")
            if "container name" in str(e) and "is already in use" in str(e):
                 console.print(f"Container name '[cyan]{CONTAINER_NAME}[/cyan]' is already in use.")
                 console.print("Try stopping or removing the existing container first ('stop' command).")
            return None
        except Exception as e:
            console.print(f"[bold red]An unexpected error occurred starting container:[/bold red] {e}")
            return None

    def stop_container(self):
        """Stops and removes the Docker container."""
        console = Console()
        container = self._find_container()
        if not container:
            console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' not found.")
            return True

        console.print(f"Stopping container '[cyan]{CONTAINER_NAME}[/cyan]' ([yellow]{container.short_id}[/yellow])...")
        stopped = False
        removed = False
        try:
            if container.status == 'running':
                container.stop()
                console.print("Container stopped.")
                stopped = True
            else:
                console.print("Container was not running.")
                stopped = True # Consider it stopped if not running

            try:
                 container.remove()
                 console.print("Container removed.")
                 removed = True
            except docker.errors.APIError as e:
                 if e.response.status_code != 404: raise e
                 else: console.print("Container was already removed.")
                 removed = True # Already gone is success

            return stopped and removed
        except Exception as e:
            console.print(f"[bold red]Error stopping/removing container:[/bold red] {e}")
            return False

    def get_status(self):
        """Gets the status of the container."""
        container = self._find_container()
        if container:
            return container.status
        return "not found"

    def run_code(self, code_string):
        """Runs a Python code string inside the running container using 'docker exec'."""
        console = Console()
        status = self.get_status()
        if status != 'running':
            console.print(f"[bold red]Error:[/bold red] Container '[cyan]{CONTAINER_NAME}[/cyan]' is not running (status: {status}).")
            console.print("Please start the container first using the 'start' command.")
            return None

        console.print(f"Executing code in container '[cyan]{CONTAINER_NAME}[/cyan]'...")
        command = ['docker', 'exec', CONTAINER_NAME, 'python3', '-c', code_string]

        try:
            result = subprocess.run(command, capture_output=True, text=True, check=False)

            console.print("[bold]--- Execution Output ---[/bold]")
            if result.stdout:
                console.print(result.stdout.strip())
            if result.stderr:
                console.print("[bold red]--- Execution Errors ---[/bold red]", file=sys.stderr)
                console.print(result.stderr.strip(), file=sys.stderr)
            console.print("[bold]------------------------[/bold]")

            if result.returncode != 0:
                console.print(f"[yellow]Code execution finished with non-zero exit code: {result.returncode}[/yellow]")
            else:
                console.print("[green]Code execution finished successfully.[/green]")
            return result.stdout # Return stdout regardless of exit code

        except FileNotFoundError:
            console.print("[bold red]Error:[/bold red] 'docker' command not found. Is Docker installed and in your PATH?")
            return None
        except Exception as e:
            console.print(f"[bold red]An error occurred running code via docker exec:[/bold red] {e}")
            return None

# --- Interactive Mode Functions ---

def print_interactive_help(console):
    """Prints help message for interactive mode."""
    console.print("\n[bold cyan]Available Commands:[/bold cyan]")
    console.print("  [green]build[/green]               Build the Docker image.")
    console.print("  [green]start[/green] [--rebuild]   Start the Docker container (optionally rebuild first).")
    console.print("  [green]stop[/green]                Stop and remove the Docker container.")
    console.print("  [green]status[/green]              Check the current status of the container.")
    console.print("  [green]run[/green] <code>          Run a Python code string inside the container.")
    console.print("                      (Example: run \"print('hello')\")")
    console.print("  [green]help[/green]                Show this help message.")
    console.print("  [green]exit[/green]                Exit the interactive manager.")
    console.print("\nExample: [yellow]start --rebuild[/yellow]")
    console.print("Example: [yellow]run \"import os; print(os.getcwd())\"[/yellow]")

def interactive_loop(manager):
    """Runs the interactive command loop."""
    console = Console()
    console.print("[bold blue]Welcome to the Interactive Benchmarking Sandbox Manager![/bold blue]")
    print_interactive_help(console)

    while True:
        try:
            if HAS_RICH:
                 raw_command = Prompt.ask("\nEnter command (\'help\' or \'exit\')")
            else:
                 raw_command = input("\nEnter command ('help' or 'exit'): ").strip()

            if not raw_command:
                continue

            try:
                command_parts = shlex.split(raw_command)
            except ValueError as e:
                console.print(f"[red]Error parsing command (check quotes?): {e}[/red]")
                continue

            if not command_parts: continue

            command = command_parts[0].lower()
            args = command_parts[1:] # Arguments passed after the command

            if command == "exit":
                # Optionally ask to stop container before exiting
                if manager.get_status() == 'running':
                     if Prompt.ask(f"Container '{CONTAINER_NAME}' is running. Stop it before exiting?", choices=["y", "n"], default="y").lower() == 'y':
                          manager.stop_container()
                break
            elif command == "help":
                print_interactive_help(console)
            elif command == "build":
                if len(args) == 0: manager.build_image()
                else: console.print("[yellow]Usage: build[/yellow]")
            elif command == "start":
                rebuild = '--rebuild' in args
                if all(a == '--rebuild' for a in args) or len(args) == 0:
                     manager.start_container(rebuild=rebuild)
                else:
                     console.print("[yellow]Usage: start [--rebuild][/yellow]")
            elif command == "stop":
                 if len(args) == 0: manager.stop_container()
                 else: console.print("[yellow]Usage: stop[/yellow]")
            elif command == "status":
                 if len(args) == 0:
                     status = manager.get_status()
                     console.print(f"Container '[cyan]{CONTAINER_NAME}[/cyan]' status: [bold {('green' if status=='running' else 'yellow')}]{status}[/bold {('green' if status=='running' else 'yellow')}]")
                 else: console.print("[yellow]Usage: status[/yellow]")
            elif command == "run":
                 if len(args) == 1:
                     code_to_run = args[0]
                     manager.run_code(code_to_run)
                 else:
                     console.print("[yellow]Usage: run \"<python_code_string>\"[/yellow]")
            else:
                console.print(f"[red]Unknown command: '{command}'. Type 'help' for options.[/red]")

        except EOFError: # Handle Ctrl+D
             console.print("\n[yellow]EOF detected. Exiting.[/yellow]")
             if manager.get_status() == 'running':
                  if Prompt.ask(f"Container '{CONTAINER_NAME}' is running. Stop it before exiting?", choices=["y", "n"], default="y").lower() == 'y':
                       manager.stop_container()
             break
        except KeyboardInterrupt: # Handle Ctrl+C
             console.print("\n[yellow]Interrupted by user. Type 'exit' to quit.[/yellow]")
             # Continue loop instead of exiting immediately
        except Exception as e:
             console.print(f"[bold red]An unexpected error occurred in the interactive loop:[/bold red] {e}")
             # Depending on severity, you might want to break or continue

    console.print("[bold blue]Exiting sandbox manager. Goodbye![/bold blue]")


# --- Main Entry Point ---
def main():
    # Initialize manager once
    try:
        manager = SandboxManager()
    except SystemExit:
        # Exit if manager initialization failed (e.g., Docker not running)
        sys.exit(1)

    # Check if running interactively (no arguments other than script name)
    if len(sys.argv) == 1:
        interactive_loop(manager)
        sys.exit(0)

    # --- Original argparse logic for non-interactive mode ---
    parser = argparse.ArgumentParser(
        description="Manage the Benchmarking Sandbox Docker container. Run without arguments for interactive mode.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Action to perform')
    subparsers.required = True # Make subcommand required if args are given

    # Build command
    parser_build = subparsers.add_parser('build', help='Build the Docker image.')
    parser_build.set_defaults(func=lambda args, mgr: mgr.build_image())

    # Start command
    parser_start = subparsers.add_parser('start', help='Start the Docker container.')
    parser_start.add_argument('--rebuild', action='store_true', help='Rebuild the image before starting.')
    parser_start.set_defaults(func=lambda args, mgr: mgr.start_container(rebuild=args.rebuild))

    # Stop command
    parser_stop = subparsers.add_parser('stop', help='Stop and remove the Docker container.')
    parser_stop.set_defaults(func=lambda args, mgr: mgr.stop_container())

    # Status command
    parser_status = subparsers.add_parser('status', help='Check if the container is running.')
    parser_status.set_defaults(func=lambda args, mgr: print(f"Container status: {mgr.get_status()}")) # Simplified output for non-interactive

    # Run command
    parser_run = subparsers.add_parser('run', help='Run a Python code string inside the container.')
    parser_run.add_argument('code', type=str, help='The Python code string to execute.')
    parser_run.set_defaults(func=lambda args, mgr: mgr.run_code(args.code))

    args = parser.parse_args()
    
    result = args.func(args, manager)
    sys.exit(0)


if __name__ == "__main__":
    main()
