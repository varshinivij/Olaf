from dataclasses import dataclass
from datetime import datetime


@dataclass
class User:
    id: str
    email: str
    name: str | None
    role: str
    status: str
    createdAt: datetime
    updatedAt: datetime
    organization: str | None
    profilePictureUrl: str | None
