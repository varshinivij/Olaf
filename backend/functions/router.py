import logging
from typing import Callable, Dict, List, Any

# Assuming Pipe is defined elsewhere and that route functions return (destination, response_generator).
from pipe import Pipe

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


class Router:
    """
    The Router directs session data through a series of transformations (pipes) 
    and then to a specified route function.

    Each route is associated with a callable that takes processed session data 
    and returns a tuple (destination, response_generator).

    The Router can:
    - Register routes and pipes (both global and route-specific).
    - Retrieve session data from a stubbed `get_session_data` method for demonstration.
    - Apply pipes to session data before routing.
    - Handle routing to multiple destinations, including chaining routes.

    Usage:
        1. Add routes with `add_route`.
        2. Add pipes if needed (global or route-specific).
        3. Call `route` with a route key and session data to process and execute.
    """

    def __init__(self):
        self.routes: Dict[str, Callable] = {}
        self.pipes: Dict[str, Pipe] = {}
        self.global_pipes: List[Pipe] = []

    def add_route(self, key: str, function: Callable) -> None:
        """
        Register a route and its corresponding function.

        Args:
            key: A unique string identifier for the route.
            function: A callable taking processed session data (dict) and returning 
                      (destination: str, response_generator: generator).
        """
        self.routes[key] = function
        logger.debug(f"Route '{key}' added.")

    def add_pipe(self, key: str, pipe: Pipe) -> None:
        """
        Associate a pipe with a specific route.

        Args:
            key: The route identifier.
            pipe: A Pipe instance to transform the session data before the route function is called.
        """
        self.pipes[key] = pipe
        logger.debug(f"Route-specific pipe added for route '{key}'.")

    def add_pipe_to_all_paths(self, pipe: Pipe) -> None:
        """
        Add a pipe that applies to all routes before route-specific pipes.

        Args:
            pipe: A Pipe instance.
        """
        self.global_pipes.append(pipe)
        logger.debug("Global pipe added for all routes.")

    def get_all_routes(self) -> Dict[str, Callable]:
        """
        Return a dictionary of all registered routes and their corresponding functions.
        """
        return self.routes

    def get_all_route_names(self) -> List[str]:
        """
        Return a list of all registered route keys.
        """
        return list(self.routes.keys())

    def route(self, key: str, session_data: Dict[str, Any]) -> Any:
        """
        Route the given session data through the appropriate pipes and to the specified route.

        Steps:
          1. Validate that the route key exists.
          2. Apply global pipes to session_data.
          3. Apply route-specific pipe if one exists for this key.
          4. Invoke the route function with the transformed session_data.
          5. Handle the returned destination:
             - If destination == 'user', return the response_generator.
             - If destination is another route key, recursively call `route` on the new key.
             - Otherwise, raise a ValueError for unknown destinations.

        Args:
            key: The route identifier.
            session_data: A dictionary containing session-related data.

        Returns:
            Any: The response from the route function or a subsequent route call.

        Raises:
            KeyError: If the specified route key doesn't exist.
            ValueError: If the destination is neither 'user' nor a known route key.
        """
        if key not in self.routes:
            logger.error(f"Route not found: {key}")
            raise KeyError(f"Route not found: {key}")

        # Apply pipes
        processed_data = self._apply_pipes(key, session_data)

        # Call the route function
        destination, response_generator = self.routes[key](processed_data)
        logger.debug(f"Route '{key}' returned destination '{destination}'.")

        # Determine next steps based on destination
        if destination == "user":
            return response_generator
        elif destination in self.routes:
            # Recursively handle routing
            return self.route(destination, processed_data)
        else:
            logger.error(f"Unknown destination: {destination}")
            raise ValueError(f"Unknown destination: {destination}")

    def route_without_pipe(self, key: str, session_data: Dict[str, Any]) -> Any:
        """
        Route a session directly to its route function without applying any pipes.
        Useful for debugging or special cases where pipe transformations are not desired.

        Args:
            key: The route identifier.
            session_id: The ID of the session to process.

        Returns:
            Any: The result from the route function directly.

        Raises:
            KeyError: If the route key doesn't exist.
            ValueError: If the destination is unknown.
        """
        if key not in self.routes:
            logger.error(f"Route not found: {key}")
            raise KeyError(f"Route not found: {key}")

        # Call the route function without any pipes
        destination, response_generator = self.routes[key](session_data)
        logger.debug(f"Route '{key}' returned destination '{destination}' without pipes.")

        if destination == "user":
            return response_generator
        elif destination in self.routes:
            return self.route(destination, session_data)
        else:
            logger.error(f"Unknown destination: {destination}")
            raise ValueError(f"Unknown destination: {destination}")

    async def _apply_pipes(self, key: str, session_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply global pipes and then the route-specific pipe (if any) to the session data asynchronously.

        Args:
            key: The route identifier for route-specific pipes.
            session_data: The initial session data dictionary.

        Returns:
            Dict[str, Any]: The transformed session data after applying all relevant pipes.
        """
        # Apply global pipes
        for pipe in self.global_pipes:
            session_data = await pipe.process(session_data)

        # Apply route-specific pipe if it exists
        if key in self.pipes:
            session_data = await self.pipes[key].process(session_data)

        return session_data