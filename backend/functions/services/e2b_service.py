from pathlib import Path
from typing import Any, Dict, List, Literal, Tuple, TypedDict

from e2b.sandbox.filesystem import FileInfo
from e2b_code_interpreter import CodeInterpreter
from e2b_code_interpreter.models import MIMEType
from google.cloud import storage
from werkzeug.datastructures.file_storage import FileStorage
from werkzeug.utils import secure_filename

from os import environ


class CodeExecutionResult(TypedDict):
    logs: Dict[Literal["stdout", "stderr"], List[str]]
    error: str | None
    results: List[Dict[MIMEType, Any]]


class TerminalExecutionResult(TypedDict):
    stdout: List[str]
    stderr: List[str]
    error: str | None


class E2BService:
    """
    Service class containing utility methods for E2B endpoints.
    """

    def __init__(self):
        # we should move this to a secret store and load from there eventually
        self.e2b_api_key = environ.get("E2B_API_KEY")
        self.e2b_template = environ.get("E2B_TEMPLATE") or "vh7kehbtf0t4xbx9ec9u"

    def create_sandbox(self) -> str:
        """
        Creates a new sandbox. Returns the sandbox id.
        """
        sandbox = CodeInterpreter(
            api_key=self.e2b_api_key,
            template=self.e2b_template,
            timeout=300,
            cwd="/home/user",
            # This is also the default folder for E2B uploads.
            # The root dir contains the Dockerfile and whatnot.
        )
        sandbox.keep_alive(5)  # keep box alive for 5 minutes
        return sandbox.id

    def get_sandbox_files(self, sandbox_id: str) -> List[FileInfo]:
        """
        Gets the files in the sandbox with the given sandbox_id.
        """
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
        files = sandbox.filesystem.list("/home/user")
        return files

    def get_is_sandbox_alive(self, sandbox_id: str) -> bool:
        """
        Gets whether the sandbox with the given sandbox_id is still alive.
        """
        try:
            sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
            return True
        except:
            return False

    def upload_formdata_files_to_sandbox(
        self, sandbox_id: str, files: List[FileStorage]
    ) -> None:
        """
        Uploads files given by an HTTP request to the sandbox with the given
        sandbox_id.

        HTTP requests should be ContentType "multipart/form-data." This can
        be done through <form> tags with enctype="multipart/form-data" or using
        FormData objects: https://developer.mozilla.org/en-US/docs/Web/API/FormData.

        Example usage: \n
        `sandbox_id = req.form["sandboxId"]` \n
        `files = req.files.getlist("file")` \n
        `upload_formdata_files_to_sandbox(sandbox_id, files)`
        """
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
        for file in files:
            if file.filename is not None:
                file.name = secure_filename(file.filename)
            sandbox.upload_file(file)  # type: ignore

    def upload_firebase_files_to_sandbox(
        self, sandbox_id: str, firebase_storage_file_paths: List[str]
    ) -> None:
        """
        Uploads files from Firebase Storage to the sandbox with the given
        sandbox_id.
        """
        bucket_ref = storage.Client().bucket("twocube-web.appspot.com")

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)

        for file_path in firebase_storage_file_paths:
            blob = bucket_ref.blob(file_path)
            file_size = blob.size
            print(f"Downloading file: {file_path}, size: {file_size} bytes")
            file_name = file_path.split("/")[-1]
            blob.download_to_filename(file_name)

            # go from file name to file object
            file_name = Path(file_name)
            with open(file_name, "rb") as f:
                remote_path = sandbox.upload_file(f)
                print(remote_path)

    def download_file_from_sandbox(
        self, sandbox_id: str, file_path: str
    ) -> Tuple[str, bytes]:
        """
        Downloads a file from the sandbox with the given sandbox_id, at a
        specified path relative to the "/home/user/" folder. Returns name of
        the file and the file bytes as a tuple.
        """
        path: Path = Path("/home/user") / file_path
        file_name = path.name

        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
        file_bytes = sandbox.download_file(path.as_posix())

        return file_name, file_bytes

    def execute_code_in_sandbox(
        self, sandbox_id: str, code: str
    ) -> CodeExecutionResult:
        """
        Returns results after executing code on the sandbox with the given
        sandbox_id.

        Returns a dict: {
            "logs": {
                "stdout": List[str] of stdout messages
                "stderr": List[str] of stderr messages
            }
            "error": None or str of E2B Error message
            "results": List[ {mimetype: str of data or None} ]
                (example.) \n
                List [{
                        "image/png": (base64 encoded PNG data as string)
                        "text/plain": "<Figure size 1000x600 with 1 Axes>"
                    },
                    ...]
            }
        """
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
        execution = sandbox.notebook.exec_cell(code)

        return {
            "logs": {"stdout": execution.logs.stdout, "stderr": execution.logs.stderr},
            "error": execution.error.traceback if execution.error else None,
            "results": [result.raw for result in execution.results],
        }

    def execute_command_in_sandbox(
        self, sandbox_id: str, command: str
    ) -> TerminalExecutionResult:
        """
        Returns results after executing a terminal command on the sandbox with
        the given sandbox_id.

        Returns a dict: {
            "stdout": List[str] of stdout messages
            "stderr": List[str] of stderr messages
            "error": None or str of E2B Error message
        }
        """
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
        execution = sandbox.notebook.exec_cell(f"!{command}")

        return {
            "stdout": execution.logs.stdout,
            "stderr": execution.logs.stderr,
            "error": execution.error.traceback if execution.error else None,
        }

    def close_sandbox(self, sandbox_id: str) -> None:
        """
        Close the sandbox with the given sandbox_id.
        """
        sandbox = CodeInterpreter.reconnect(sandbox_id, api_key=self.e2b_api_key)
        sandbox.close()
