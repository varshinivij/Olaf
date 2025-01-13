from typing import List
from models.chat_message import ChatMessage  # or wherever ChatMessage is defined
from agents.l3_coder_agent import CoderAgent  # adjust import path to your setup

def l3_coderagent_route(plan: str):
    """
    Route function for the L3 Coder Agent. 
    This route is used to generate text responses (and call tools as needed).
    """
    coder_agent = CoderAgent("python", plan)  # type: ignore
    return coder_agent.generate_response()