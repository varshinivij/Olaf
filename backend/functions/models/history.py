from typing import List
import json

from functions.models.chat_message import ChatMessage, ChatMessageRole, ChatMessageType


class History:
    def __init__(self, system: str | List[ChatMessage]):
        if isinstance(system, str):
            self.history = json.loads(system)[0]
        else:
            self.history = system

    def log(self, role: ChatMessageRole, content: str, type: ChatMessageType):
        entry: ChatMessage = {"role": role, "content": content, "type": type}
        self.history.append(entry)

    def get_history(self):
        return self.history

    def most_recent_entry(self):
        return self.history[-1]

    def remove_system_messages(self):
        self.history = [entry for entry in self.history if entry["role"] != "system"]

    def upsert(self, prompt: str, type: ChatMessageType):
        self.history.insert(0, {"role": "system", "content": prompt, "type": type})
