"""Shared Pydantic request/response models for the API."""

from pydantic import BaseModel, EmailStr, Field, field_validator

import opencloning_linkml.datamodel.models as opencloning_models
from opencloning_db.models import SequenceType, Sequence


# --- Auth (OAuth2 password + JWT) ---
class Token(BaseModel):
    access_token: str
    token_type: str = 'bearer'


class UserPublic(BaseModel):
    id: int
    email: str
    display_name: str | None
    is_instance_admin: bool


class WorkspaceRef(BaseModel):
    id: int
    name: str
    role: str


class WorkspaceCreate(BaseModel):
    name: str = Field(min_length=1)


class WorkspaceRename(BaseModel):
    name: str = Field(min_length=1)


class RegisterBody(BaseModel):
    email: EmailStr
    password: str = Field(min_length=1)
    display_name: str | None = None


# --- Sequence sample ---
class SequenceSampleCreate(BaseModel):
    uid: str
    sequence_id: int


class SequenceSampleUpdate(BaseModel):
    sequence_id: int


class SequenceSampleRead(BaseModel):
    id: int
    uid: str
    sequence_id: int


class SequenceSampleCreated(BaseModel):
    id: int
    uid: str


# --- Tags ---
class TagCreate(BaseModel):
    name: str = Field(min_length=1)

    @field_validator('name', mode='before')
    @classmethod
    def strip_tag_name(cls, v: object) -> object:
        # We do it before to strip before counting the length of the string
        if isinstance(v, str):
            return v.strip()
        return v


class TagRead(BaseModel):
    id: int
    name: str


class EntityTagAttach(BaseModel):
    tag_id: int


# --- Entity refs ---
class InputEntityRef(BaseModel):
    id: int
    type: str
    name: str | None


class SequencingFileRef(BaseModel):
    id: int
    original_name: str


# --- Common responses ---
class IdResponse(BaseModel):
    id: int


class RemovedResponse(BaseModel):
    removed: int


class DeletedResponse(BaseModel):
    deleted: int


# --- Cloning strategy ---
class CloningStrategyIdMapping(BaseModel):
    localId: int
    databaseId: int


class CloningStrategyResponse(BaseModel):
    id: int
    mappings: list[CloningStrategyIdMapping]


# --- Sequence / primer refs ---
class SequenceRef(BaseModel):
    id: int
    name: str | None
    sequence_type: SequenceType
    tags: list[TagRead] = []
    sample_uids: list[str] = []
    seguid: str | None = None


class SequenceUpdate(BaseModel):
    name: str | None = None
    sequence_type: SequenceType | None = None


class PrimerUpdate(BaseModel):
    name: str = Field(min_length=2)


class PrimerRef(BaseModel):
    id: int
    name: str | None
    sequence: str
    uid: str | None = None
    tags: list[TagRead] = []


class SequenceSampleWithSequence(BaseModel):
    id: int
    uid: str
    sequence_id: int
    sequence: SequenceRef


# --- Line ---
class SequenceInLineRef(BaseModel):
    """Sequence in a line, including the SequenceInLine instance id."""

    id: int
    sequence_id: int
    name: str | None
    sequence_type: SequenceType
    tags: list[TagRead]


class LineRef(BaseModel):
    id: int
    uid: str
    sequences_in_line: list[SequenceInLineRef]
    parent_ids: list[int]
    tags: list[TagRead]


class LineCreate(BaseModel):
    uid: str
    allele_ids: list[int] = []
    plasmid_ids: list[int] = []
    parent_ids: list[int] = []


class LineUpdateLinks(BaseModel):
    allele_ids: list[int] | None = None
    plasmid_ids: list[int] | None = None
    parent_ids: list[int] | None = None


def sequence_ref(sequence: Sequence) -> SequenceRef:
    return SequenceRef(
        id=sequence.id,
        name=sequence.name,
        sequence_type=sequence.sequence_type,
        tags=[TagRead(id=t.id, name=t.name) for t in sequence.tags],
        sample_uids=sequence.sample_uids,
        seguid=sequence.seguid,
    )


class SequenceSearchResult(BaseModel):
    sequence_ref: SequenceRef
    sequence: opencloning_models.TextFileSequence
    shift: int = 0
    reverse_complemented: bool = False
