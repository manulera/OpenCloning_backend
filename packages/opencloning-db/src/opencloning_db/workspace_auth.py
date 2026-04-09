"""Workspace membership checks for API routes."""

from sqlalchemy import select
from sqlalchemy.orm import Session

from fastapi import HTTPException, status

from opencloning_db.models import WorkspaceMembership, WorkspaceRole

_ROLE_ORDER: dict[WorkspaceRole, int] = {
    WorkspaceRole.viewer: 0,
    WorkspaceRole.editor: 1,
    WorkspaceRole.owner: 2,
}


def has_at_least(role: WorkspaceRole, minimum: WorkspaceRole) -> bool:
    return _ROLE_ORDER[role] >= _ROLE_ORDER[minimum]


def assert_workspace_access(
    session: Session,
    user_id: int,
    workspace_id: int,
    min_role: WorkspaceRole,
) -> WorkspaceMembership:
    """
    Require membership in workspace with at least ``min_role``.
    Raises 403 otherwise.
    """
    membership = session.scalar(
        select(WorkspaceMembership).where(
            WorkspaceMembership.user_id == user_id,
            WorkspaceMembership.workspace_id == workspace_id,
        )
    )
    if membership is None:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail='Not allowed for this workspace')
    if not has_at_least(membership.role, min_role):
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail='Not allowed for this workspace')
    return membership
