from typing import Literal, TypedDict


type ChatMessageType = Literal[
    "text", "code", "executedCode", "plan", "error", "image", "result", "hidden"
]
type ChatMessageRole = Literal["assistant", "user"]


class ChatMessage(TypedDict):
    """
    ChatMessages are JSON objects stored as dictionaries.
    """

    type: ChatMessageType
    role: ChatMessageRole
    content: str
