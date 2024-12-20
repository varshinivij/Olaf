import logging
from typing import List, Callable, Dict, Any

logger = logging.getLogger(__name__)

class Pipe:
    """
    A Pipe represents a pipeline of asynchronous operations to be performed on session data.

    Each operation in `operations` is an async callable that takes a session dictionary
    and returns a modified session dictionary. If any operation fails, the entire pipeline fails.

    Usage:
        1. Instantiate a Pipe with a name and a list of async operations.
        2. Call `await pipe.process(session_data)` to apply all operations in order.
    """

    def __init__(self, name: str, operations: List[Callable[[Dict[str, Any]], Any]]):
        """
        Initialize the Pipe.

        Args:
            name: A descriptive name for this pipeline of operations (e.g., "PreprocessingPipe").
            operations: A list of async functions. Each function should accept a session dict and 
                        return a (possibly modified) session dict.
        """
        self.name = name
        self.operations = operations if operations is not None else []

    def process(self, session: Dict[str, Any]) -> Dict[str, Any]:
        """
        Asynchronously process the session data through each operation in sequence.

        If any operation raises an exception, the process is halted and the exception is re-raised.

        Args:
            session: A dictionary containing session-related data to transform.

        Returns:
            A dictionary representing the transformed session data after all operations have been applied.

        Raises:
            Exception: If any operation fails, the original exception is logged and re-raised.
        """
        processed_session = session.copy()
        for idx, operation in enumerate(self.operations):
            try:
                processed_session = operation(processed_session)
            except Exception as e:
                logger.exception(f"Error in pipe '{self.name}' at operation index {idx}: {str(e)}")
                raise
        return processed_session