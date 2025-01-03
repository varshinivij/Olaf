from abc import ABC, abstractmethod
from typing import Any, Dict, List, Callable, Generator, Optional, Tuple


class AbstractAgent(ABC):
    """
    A standard abstract base class for an agent in our multi-agent AI system.

    This class defines the essential interface and behavior that all
    multi-agent controllers should implement. It provides a blueprint
    for:
      - Managing system and agent-specific prompts.
      - Handling conversation history.
      - Integrating with function/tool APIs.
      - Generating responses from an underlying LLM or reasoning engine.
    """

    def __init__(
        self,
        system_prompt: str,
        history: Any,
        functions: List[Dict[str, Any]] = [],
    ):
        """
        Initialize the multi-agent system.

        Args:
            system_prompt (str): The initial system prompt guiding the entire system.
            history (Any): An object representing the conversation or event history.
                           Could be a list of messages, a specialized history object, etc.
            functions (List[Dict[str, Any]], optional): A list of function specifications
                                                       (tools) available to agents.
        """
        self.system_prompt = system_prompt
        self.history = history
        self.functions = functions
        self.function_map = self._build_function_map()

    @abstractmethod
    def _build_function_map(self) -> Dict[str, Callable]:
        """
        Build a mapping from function/tool names to their callables.

        Returns:
            Dict[str, Callable]: A dictionary mapping function/tool names (str) to
                                 corresponding callables that implement the tool logic.
        """
        pass

    @abstractmethod
    def generate_response(self) -> Tuple[str, Generator[Dict[str, Any], None, None]]:
        """
        Generate a response to the user request.

        Returns:
            (destination, generator):
                destination: A string indicating where the response should be directed (e.g., "user").
                generator: A generator that yields dictionaries representing messages.
        """
        pass

    @abstractmethod
    def handle_functions(self, function_name: str, arguments: Dict[str, Any]) -> Any:
        """
        Call the specified function/tool with the provided arguments.

        Implementations should:
          - Look up the function in `self.function_map`.
          - Execute it safely.
          - Handle exceptions and return results in a standardized format.

        Args:
            function_name (str): The name of the function/tool to call.
            arguments (Dict[str, Any]): The arguments to pass to the function.

        Returns:
            Any: The result from the function call (could be structured data).
        """
        pass

    @abstractmethod
    def store_interaction(self, role: str, content: str) -> None:
        """
        Store the interaction (e.g., a user query or agent response) in the system's history.

        Args:
            role (str): The role of the message sender (e.g., "user", "assistant", "system").
            content (str): The content of the message.
        """
        pass

    def add_function_spec(self, function_spec: Dict[str, Any]) -> None:
        """
        Add a new function/tool specification to the system.

        Args:
            function_spec (Dict[str, Any]): A dictionary describing the function (name, description, parameters).
        """
        self.functions.append(function_spec)
        # Optionally rebuild the function map if needed
        self.function_map = self._build_function_map()
