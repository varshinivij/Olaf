from dataclasses import dataclass
from datetime import datetime
from typing import Literal


type ProjectLanguage = Literal["Python"]
type ProjectModel = Literal["GPT-4o"]


@dataclass
class Project:
    id: str
    name: str
    language: ProjectLanguage
    model: ProjectModel
    updatedAt: datetime
