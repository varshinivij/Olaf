from typing import List
from datastructures.history import History
from models.chat_message import ChatMessage  # or wherever ChatMessage is defined
from agents.l3_master_agent import MasterAgent  # adjust import path to your setup

def l3_masteragent_route(session_history: List[ChatMessage]):
    """
    Route function for the L3 Master Agent. 
    This route is used to generate text responses (and call tools as needed).
    """
    history = History(session_history)
    master_agent = MasterAgent("python", history)  # type: ignore
    return master_agent.generate_response()