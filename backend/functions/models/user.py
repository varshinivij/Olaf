from datetime import datetime
from typing import TypedDict


class User(TypedDict):
    """
    Users are JSON objects stored as dictionaries.
    """

    id: str
    email: str
    name: str | None
    role: str
    status: str
    createdAt: datetime
    updatedAt: datetime
    organization: str | None
    profilePictureUrl: str | None
