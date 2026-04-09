"""Primer endpoints."""

from typing import Annotated

import opencloning_linkml.datamodel.models as opencloning_models
from fastapi import APIRouter, Depends, Query
from sqlalchemy import and_, select
from sqlalchemy.orm import selectinload

from opencloning_db.apimodels import IdResponse, PrimerRef, SequenceRef, TagRead, sequence_ref, PrimerUpdate
from fastapi_pagination import Page
from fastapi_pagination.ext.sqlalchemy import paginate
from opencloning_db.models import InputEntity, Primer, Sequence, Source, SourceInput, Tag, WorkspaceRole
from pydantic import create_model
from opencloning_db.workspace_deps import (
    WorkspaceContext,
    get_editor_workspace_ctx,
    get_primer_in_workspace_for_user,
    get_viewer_workspace_ctx,
)

router = APIRouter(tags=['primers'])


@router.get('/primers', response_model=Page[PrimerRef])
def get_primers(
    ctx: Annotated[WorkspaceContext, Depends(get_viewer_workspace_ctx)],
    tags: list[int] = Query(description='Filter primers by tag IDs', default_factory=list),
    name: str | None = Query(description='Filter primers by name (case-insensitive substring match)', default=None),
    uid: str | None = Query(
        description='Filter primers by sample uid (case-insensitive substring match)', default=None
    ),
    has_uid: bool = Query(description='Filter primers by whether they have a uid', default=False),
):
    current_user, session, workspace_id = ctx

    query = select(Primer).where(Primer.workspace_id == workspace_id).options(selectinload(InputEntity.tags))
    if tags:
        query = query.where(InputEntity.tags.any(and_(Tag.id.in_(tags), Tag.workspace_id == workspace_id)))
    if name is not None:
        query = query.where(Primer.name.ilike(f"%{name}%"))
    if uid is not None:
        query = query.where(Primer.uid.ilike(f"%{uid}%"))
    if has_uid is True:
        query = query.where(Primer.uid.isnot(None))
    return paginate(session, query)


@router.post('/primer', response_model=IdResponse)
def post_primer(
    ctx: Annotated[WorkspaceContext, Depends(get_editor_workspace_ctx)],
    primer: opencloning_models.Primer,
):
    """Submit a standalone primer (unlinked to any cloning strategy)."""
    current_user, session, workspace_id = ctx

    db_primer = Primer.from_pydantic(primer, workspace_id)
    session.add(db_primer)
    session.commit()
    session.refresh(db_primer)
    return IdResponse(id=db_primer.id)


@router.get('/primer/{primer_id}', response_model=PrimerRef)
def get_primer(
    primer_id: int,
    ctx: Annotated[WorkspaceContext, Depends(get_viewer_workspace_ctx)],
):
    current_user, session, workspace_id = ctx
    primer = get_primer_in_workspace_for_user(session, current_user, workspace_id, primer_id, WorkspaceRole.viewer)

    return PrimerRef(
        id=primer.id,
        name=primer.name,
        sequence=primer.sequence,
        uid=primer.uid,
        tags=[TagRead(id=t.id, name=t.name) for t in primer.tags],
    )


@router.get(
    '/primer/{primer_id}/sequences',
    response_model=create_model('PrimerSequences', templates=list[SequenceRef], products=list[SequenceRef]),
)
def get_primer_sequences(
    primer_id: int,
    ctx: Annotated[WorkspaceContext, Depends(get_viewer_workspace_ctx)],
):
    """Get sequences linked to a primer."""
    current_user, session, workspace_id = ctx
    # Check that the user has access to the primer and the primer exists
    get_primer_in_workspace_for_user(session, current_user, workspace_id, primer_id, WorkspaceRole.viewer)

    # 1) All source IDs where this primer was an input
    primer_source_ids = select(SourceInput.source_id).where(SourceInput.input_entity_id == primer_id)
    # 2) All sequences that were inputs to those same sources
    template_sequences = (
        select(Sequence)
        .join(
            SourceInput,
            SourceInput.input_entity_id == Sequence.id,
        )
        .where(
            SourceInput.source_id.in_(primer_source_ids),
            Sequence.workspace_id == Primer.workspace_id,  # Safety check TODO: Maybe remove?
        )
        .distinct()
    )
    # 3) All sequences that were outputs of those same sources
    product_sequences = (
        select(Sequence)
        .join(Source)
        .where(
            Source.id.in_(primer_source_ids),
            Sequence.workspace_id == Primer.workspace_id,  # Safety check TODO: Maybe remove?
        )
    )
    templates = session.execute(template_sequences).scalars().all()
    products = session.execute(product_sequences).scalars().all()

    return {
        'templates': [sequence_ref(s) for s in templates],
        'products': [sequence_ref(s) for s in products],
    }


@router.patch('/primer/{primer_id}', response_model=PrimerRef)
def patch_primer(
    primer_id: int,
    body: PrimerUpdate,
    ctx: Annotated[WorkspaceContext, Depends(get_editor_workspace_ctx)],
):
    """
    Update primer name.

    Note: the primer "type" is fixed by polymorphic identity and cannot be changed.
    """
    current_user, session, workspace_id = ctx
    primer = get_primer_in_workspace_for_user(session, current_user, workspace_id, primer_id, WorkspaceRole.editor)
    primer.name = body.name

    session.commit()
    session.refresh(primer)

    return PrimerRef(
        id=primer.id,
        name=primer.name,
        sequence=primer.sequence,
        uid=primer.uid,
        tags=[TagRead(id=t.id, name=t.name) for t in primer.tags],
    )
