"""Workspace listing endpoints for the current user."""

from typing import Annotated

from fastapi import APIRouter, Depends
from sqlalchemy import select
from sqlalchemy.orm import Session

from apimodels import WorkspaceCreate, WorkspaceRef, WorkspaceRename
from deps import get_current_user, get_db
from models import User, Workspace, WorkspaceMembership, WorkspaceRole
from workspace_auth import assert_workspace_access

router = APIRouter(tags=['workspaces'])


@router.post('/workspaces', response_model=WorkspaceRef)
def create_workspace(
    body: WorkspaceCreate,
    current_user: Annotated[User, Depends(get_current_user)],
    session: Annotated[Session, Depends(get_db)],
) -> WorkspaceRef:
    workspace = Workspace(name=body.name)
    session.add(workspace)
    session.flush()
    session.add(
        WorkspaceMembership(
            user_id=current_user.id,
            workspace_id=workspace.id,
            role=WorkspaceRole.owner,
        )
    )
    session.commit()
    return WorkspaceRef(
        id=workspace.id,
        name=workspace.name,
        role=WorkspaceRole.owner.value,
    )


@router.get('/workspaces', response_model=list[WorkspaceRef])
def list_workspaces(
    current_user: Annotated[User, Depends(get_current_user)],
    session: Annotated[Session, Depends(get_db)],
) -> list[WorkspaceRef]:
    memberships = session.scalars(
        select(WorkspaceMembership)
        .join(Workspace, Workspace.id == WorkspaceMembership.workspace_id)
        .where(WorkspaceMembership.user_id == current_user.id)
        .order_by(Workspace.id.asc())
    ).all()
    return [
        WorkspaceRef(
            id=membership.workspace.id,
            name=membership.workspace.name,
            role=membership.role.value,
        )
        for membership in memberships
    ]


@router.get('/workspaces/{workspace_id}', response_model=WorkspaceRef)
def get_workspace(
    workspace_id: int,
    current_user: Annotated[User, Depends(get_current_user)],
    session: Annotated[Session, Depends(get_db)],
) -> WorkspaceRef:
    membership = assert_workspace_access(
        session,
        current_user.id,
        workspace_id,
        WorkspaceRole.viewer,
    )
    workspace = session.get(Workspace, workspace_id)
    return WorkspaceRef(
        id=workspace.id,
        name=workspace.name,
        role=membership.role.value,
    )


@router.patch('/workspaces/{workspace_id}', response_model=WorkspaceRef)
def rename_workspace(
    workspace_id: int,
    body: WorkspaceRename,
    current_user: Annotated[User, Depends(get_current_user)],
    session: Annotated[Session, Depends(get_db)],
) -> WorkspaceRef:
    membership = assert_workspace_access(
        session,
        current_user.id,
        workspace_id,
        WorkspaceRole.owner,
    )
    workspace = session.get(Workspace, workspace_id)
    workspace.name = body.name
    session.commit()
    return WorkspaceRef(
        id=workspace.id,
        name=workspace.name,
        role=membership.role.value,
    )
