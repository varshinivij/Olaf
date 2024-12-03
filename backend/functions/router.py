from urllib import response
from pipe import Pipe
import firebase_admin
from firebase_admin import credentials, firestore, auth
from functools import wraps
from typing import Callable, Dict, List, Any, Optional
import time
from datetime import datetime

class Router:
    def __init__(self):
        self.routes: Dict[str, Callable] = {}
        self.pipes: Dict[str, Pipe] = {}
        self.global_pipes: List[Pipe] = []

    def add_route(self, key: str, function: Callable) -> None:
        self.routes[key] = function

    def add_pipe(self, key: str, pipe: Pipe) -> None:
        self.pipes[key] = pipe

    def add_pipe_to_all_paths(self, pipe: Pipe) -> None:
        self.global_pipes.append(pipe)

    def get_all_routes(self) -> Dict[str, Callable]:
        return self.routes

    def get_all_route_names(self) -> List[str]:
        return list(self.routes.keys())

    def get_session_data(self, session_id: str) -> dict:
        """
        Return dummy session data with the specified fields, including a history of LLM-like responses.
        
        Args:
            session_id (str): The ID of the session to fetch
            
        Returns:
            dict: Dummy session data including context, history, and other fields.
        """
        # Dummy session data structure
        return {
            "context": "",
            "files": [],
            "history": [
                {"role": "user", "content": "Hello, what can you do?"},
                {"role": "assistant", "content": "I can help you with various tasks such as answering questions, writing code, and more."},
                {"role": "user", "content": "Can you write a function to calculate average in python?"},
            ],
            "name": f"Session {session_id}",
            "sandboxId": None,
            "id":"0f3dKnxE0AJTdXqMlive",
            "userId": "ICkQVNZhmhZkfHt00s4gCOU3BQu2"
        }

    def route_session(self, key: str, session) -> Any:
        """
        Route a session through pipes and to its destination function.
        
        Args:
            key (str): Route identifier
            session_id (str): ID of the session to process
            
        Returns:
            Any: Result from the route function
        """
        if key not in self.routes:
            raise KeyError(f"Route not found: {key}")
                
        processed_session = session.copy()
        # Apply global pipes
        for pipe in self.global_pipes:
            processed_session = pipe.process(processed_session)
        # Apply route-specific pipe if it exists
        if key in self.pipes:
            processed_session = self.pipes[key].process(processed_session)
        # Get the destination and response generator
        destination, response_generator = self.routes[key](processed_session)
        print("destination: ", destination)
        # Handle routing based on the destination
        if destination == "user":
            return response_generator
        elif destination in self.routes:
            # Route to another agent
            return self.route_session(destination, processed_session)
        else:
            raise ValueError(f"Unknown destination: {destination}")

    async def route_without_pipe(self, key: str, session_id: str) -> Any:
        """
        Route a session directly to its function without applying pipes.
        
        Args:
            key (str): Route identifier
            session_id (str): ID of the session to process
            
        Returns:
            Any: Result from the route function
        """
        if key not in self.routes:
            raise KeyError(f"Route not found: {key}")
            
        session = await self.get_session_data(session_id)
        return await self.routes[key](session)