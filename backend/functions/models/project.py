from datetime import datetime
from typing import Literal, TypedDict


type ProjectLanguage = Literal["Python"]
type ProjectModel = Literal["GPT-4o"]


class Project(TypedDict):
    """
    Projects are JSON objects stored as dictionaries.
    """

    id: str
    name: str
    language: ProjectLanguage
    model: ProjectModel
    updatedAt: datetime
