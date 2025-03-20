from typing import List

from agents.codemaster_agent import CodeMasterAgent
from datastructures.history import History
from models.chat_message import ChatMessage


def masteragent_route(session_history: List[ChatMessage]):
    """
    Master route for CodeMasterAgent. This route is used to generate text and
    code responses.
    """
    history = History(session_history)
    master_agent = CodeMasterAgent("python", history)  # type: ignore
    return master_agent.generate_response()
