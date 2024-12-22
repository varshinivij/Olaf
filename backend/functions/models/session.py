from typing import List, TypedDict
from chat_message import ChatMessage


class Session(TypedDict):
    """
    Sessions are JSON objects stored as dictionaries.
    """

    id: str
    userId: str
    projectId: str
    name: str
    context: str
    history: List[ChatMessage]
    sandboxId: str | None
