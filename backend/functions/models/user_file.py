from datetime import datetime
from typing import TypedDict


class UserFile(TypedDict):
    """
    UserFiles are JSON objects stored as dictionaries.
    """

    id: str
    name: str
    path: str
    size: int
    extension: str
    isFolder: bool
    uploadedOn: datetime
    storageLink: str
