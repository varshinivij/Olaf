from typing import List
import json

from functions.models.chat_message import ChatMessage, ChatMessageRole, ChatMessageType


class History:
    """
    Utility class for handling chat history. Wraps a list of ChatMessages.
    """

    def __init__(self, system: str | List[ChatMessage]):
        if isinstance(system, str):
            self.history: List[ChatMessage] = json.loads(system)[0]
        else:
            self.history: List[ChatMessage] = system

    def log(self, role: ChatMessageRole, content: str, type: ChatMessageType) -> None:
        entry: ChatMessage = {"role": role, "content": content, "type": type}
        self.history.append(entry)

    def get_history(self) -> List[ChatMessage]:
        return self.history

    def most_recent_entry(self) -> ChatMessage:
        return self.history[-1]

    def remove_system_messages(self) -> None:
        self.history = [entry for entry in self.history if entry["role"] != "system"]

    def upsert(self, prompt: str, type: ChatMessageType) -> None:
        self.history.insert(0, {"role": "system", "content": prompt, "type": type})
