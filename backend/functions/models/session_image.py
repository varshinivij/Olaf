from datetime import datetime
from typing import TypedDict


class SessionImage(TypedDict):
    """
    SessionImages are JSON objects stored as dictionaries.
    """

    id: str
    name: str
    path: str
    size: int
    storageLink: str
    uploadedOn: datetime
    user_id: str
