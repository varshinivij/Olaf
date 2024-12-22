from dataclasses import dataclass
from typing import Literal


@dataclass
class ChatMessage:
    type: Literal[
        "text", "code", "executedCode", "plan", "error", "image", "result", "hidden"
    ]
    role: Literal["assistant", "user"]
    content: str
