"""Shared test utilities for workspace resource fixtures."""

import sys

import config as config_module
from config import set_config
from auth.security import get_password_hash

_JWT_SECRET = 'test-jwt-secret-not-for-production'

_STANDARD_EMAILS = {
    'owner_w1': 'owner-w1@test.com',
    'owner_w2': 'owner-w2@test.com',
    'viewer_w1': 'viewer-w1@test.com',
    'owner_both': 'owner-both@test.com',
    'owner_w1_viewer_w2': 'owner-w1-viewer-w2@test.com',
}

_STANDARD_PASSWORDS = {
    'owner_w1': 'pw-owner-w1',
    'owner_w2': 'pw-owner-w2',
    'viewer_w1': 'pw-viewer-w1',
    'owner_both': 'pw-owner-both',
    'owner_w1_viewer_w2': 'pw-owner-w1-viewer-w2',
}


def bearer_headers(token: str) -> dict[str, str]:
    """Return ``Authorization: Bearer`` header only."""
    return {'Authorization': f"Bearer {token}"}


def workspace_headers(
    token: str,
    workspace_id: int | str,
    extra: dict[str, str] | None = None,
) -> dict[str, str]:
    """Bearer plus ``X-Workspace-Id`` (coerced to str).

    *extra* merges additional headers (e.g. ``Content-Type`` for JSON bodies).
    """
    h = {
        **bearer_headers(token),
        'X-Workspace-Id': str(workspace_id),
    }
    if extra:
        h.update(extra)
    return h


def post_sequencing_file_upload(
    client,
    seq_id: int,
    token: str,
    workspace_id: int | str,
    filename: str,
    body: bytes,
):
    """POST multipart upload to ``/sequence/{id}/sequencing_files``."""
    return client.post(
        f"/sequence/{seq_id}/sequencing_files",
        headers=workspace_headers(token, workspace_id),
        files={'files': (filename, body, 'application/octet-stream')},
    )


def assert_get_missing_workspace_header_422(client, path: str, token: str) -> None:
    """GET with bearer only (no X-Workspace-Id) returns 422."""
    r = client.get(path, headers=bearer_headers(token))
    assert r.status_code == 422, r.text


def assert_get_invalid_workspace_id_422(
    client,
    path: str,
    token: str,
    invalid: str = 'bad',
) -> None:
    """GET with non-coercible X-Workspace-Id returns 422."""
    r = client.get(
        path,
        headers={**bearer_headers(token), 'X-Workspace-Id': invalid},
    )
    assert r.status_code == 422, r.text


def assert_get_non_member_workspace_403(
    client,
    path: str,
    token: str,
    workspace_id: int | str,
) -> None:
    """GET with bearer and workspace id user cannot access returns 403."""
    r = client.get(path, headers=workspace_headers(token, workspace_id))
    assert r.status_code == 403, r.text


def assert_get_unauthenticated_401(
    client,
    path: str,
    workspace_id: int | str,
) -> None:
    """GET with X-Workspace-Id only (no bearer) returns 401."""
    r = client.get(path, headers={'X-Workspace-Id': str(workspace_id)})
    assert r.status_code == 401, r.text


def assert_post_unauthenticated_401(
    client,
    path: str,
    workspace_id: int | str,
    *,
    json,
) -> None:
    """POST with X-Workspace-Id only (no bearer) returns 401."""
    r = client.post(
        path,
        headers={'X-Workspace-Id': str(workspace_id)},
        json=json,
    )
    assert r.status_code == 401, r.text


def assert_patch_unauthenticated_401(
    client,
    path: str,
    workspace_id: int | str,
    *,
    json,
) -> None:
    """PATCH with X-Workspace-Id only (no bearer) returns 401."""
    r = client.patch(
        path,
        headers={'X-Workspace-Id': str(workspace_id)},
        json=json,
    )
    assert r.status_code == 401, r.text


def make_app_client(
    tmp_path,
    monkeypatch,
    extra_module: str,
    *,
    sequence_files_dir: str | None = None,
    sequencing_files_dir: str | None = None,
):
    """Bootstrap a fresh SQLite DB and return ``(engine, TestClient)``.

    Purges ``sys.modules`` for the standard set plus *extra_module* (the
    router under test) so each fixture starts from a clean import state.

    Optional *sequence_files_dir* and *sequencing_files_dir* override config
    paths (e.g. under ``tmp_path``) so tests can control on-disk GenBank and
    sequencing uploads without touching the repo default directories.
    """
    monkeypatch.setenv('OPENCLONING_JWT_SECRET', _JWT_SECRET)
    db_path = tmp_path / 'test.db'
    cfg_kwargs: dict = {
        'database_url': f"sqlite:///{db_path}",
        'jwt_secret': _JWT_SECRET,
    }
    if sequence_files_dir is not None:
        cfg_kwargs['sequence_files_dir'] = sequence_files_dir
    if sequencing_files_dir is not None:
        cfg_kwargs['sequencing_files_dir'] = sequencing_files_dir
    set_config(config_module.Config(**cfg_kwargs))

    for name in list(sys.modules.keys()):
        if name in ('api', 'deps', 'routers.auth', 'db', extra_module):
            del sys.modules[name]

    import db as db_module
    from fastapi.testclient import TestClient
    from models import Base

    engine = db_module.get_engine(config_module.get_config())
    Base.metadata.create_all(engine)

    from api import app

    return engine, TestClient(app)


def seed_standard_users(session) -> dict:
    """Seed the standard users / workspaces topology and return a context dict.

    Membership topology:
        owner_w1            — owner of W1 only
        owner_w2            — owner of W2 only
        viewer_w1           — viewer of W1 only
        owner_both          — owner of both W1 and W2
        owner_w1_viewer_w2  — owner of W1, viewer of W2

    W3 exists but has no standard memberships (for forbidden-access tests).

    The session is flushed but **not** committed; callers should add their
    own resource rows and commit once.
    """
    from models import User, Workspace, WorkspaceMembership, WorkspaceRole

    owner_w1 = User(
        email=_STANDARD_EMAILS['owner_w1'],
        display_name='Owner W1',
        password_hash=get_password_hash(_STANDARD_PASSWORDS['owner_w1']),
    )
    owner_w2 = User(
        email=_STANDARD_EMAILS['owner_w2'],
        display_name='Owner W2',
        password_hash=get_password_hash(_STANDARD_PASSWORDS['owner_w2']),
    )
    viewer_w1 = User(
        email=_STANDARD_EMAILS['viewer_w1'],
        display_name='Viewer W1',
        password_hash=get_password_hash(_STANDARD_PASSWORDS['viewer_w1']),
    )
    owner_both = User(
        email=_STANDARD_EMAILS['owner_both'],
        display_name='Owner Both',
        password_hash=get_password_hash(_STANDARD_PASSWORDS['owner_both']),
    )
    owner_w1_viewer_w2 = User(
        email=_STANDARD_EMAILS['owner_w1_viewer_w2'],
        display_name='Owner W1 Viewer W2',
        password_hash=get_password_hash(_STANDARD_PASSWORDS['owner_w1_viewer_w2']),
    )
    w1 = Workspace(name='Workspace One')
    w2 = Workspace(name='Workspace Two')
    w3 = Workspace(name='Workspace Three')
    session.add_all([owner_w1, owner_w2, viewer_w1, owner_both, owner_w1_viewer_w2, w1, w2, w3])
    session.flush()

    session.add_all(
        [
            WorkspaceMembership(
                user_id=owner_w1.id,
                workspace_id=w1.id,
                role=WorkspaceRole.owner,
            ),
            WorkspaceMembership(
                user_id=owner_w2.id,
                workspace_id=w2.id,
                role=WorkspaceRole.owner,
            ),
            WorkspaceMembership(
                user_id=viewer_w1.id,
                workspace_id=w1.id,
                role=WorkspaceRole.viewer,
            ),
            WorkspaceMembership(
                user_id=owner_both.id,
                workspace_id=w1.id,
                role=WorkspaceRole.owner,
            ),
            WorkspaceMembership(
                user_id=owner_both.id,
                workspace_id=w2.id,
                role=WorkspaceRole.owner,
            ),
            WorkspaceMembership(
                user_id=owner_w1_viewer_w2.id,
                workspace_id=w1.id,
                role=WorkspaceRole.owner,
            ),
            WorkspaceMembership(
                user_id=owner_w1_viewer_w2.id,
                workspace_id=w2.id,
                role=WorkspaceRole.viewer,
            ),
        ]
    )

    return {
        'w1': w1.id,
        'w2': w2.id,
        'w3': w3.id,
        'owner_w1_id': owner_w1.id,
        'owner_w1_email': _STANDARD_EMAILS['owner_w1'],
        'owner_w1_pw': _STANDARD_PASSWORDS['owner_w1'],
        'owner_w2_id': owner_w2.id,
        'owner_w2_email': _STANDARD_EMAILS['owner_w2'],
        'owner_w2_pw': _STANDARD_PASSWORDS['owner_w2'],
        'viewer_w1_id': viewer_w1.id,
        'viewer_w1_email': _STANDARD_EMAILS['viewer_w1'],
        'viewer_w1_pw': _STANDARD_PASSWORDS['viewer_w1'],
        'owner_both_id': owner_both.id,
        'owner_both_email': _STANDARD_EMAILS['owner_both'],
        'owner_both_pw': _STANDARD_PASSWORDS['owner_both'],
        'owner_w1_viewer_w2_id': owner_w1_viewer_w2.id,
        'owner_w1_viewer_w2_email': _STANDARD_EMAILS['owner_w1_viewer_w2'],
        'owner_w1_viewer_w2_pw': _STANDARD_PASSWORDS['owner_w1_viewer_w2'],
    }


def fetch_token(client, email: str, password: str) -> str:
    """Obtain a JWT access token for the given credentials."""
    r = client.post(
        '/auth/token',
        data={'username': email, 'password': password},
    )
    assert r.status_code == 200, r.text
    return r.json()['access_token']


def attach_standard_tokens(ctx: dict, client) -> dict:
    """Fetch tokens for all standard seeded users and attach them to *ctx*.

    Also sets ``ctx["client"]``. Returns *ctx* for convenience.
    """
    ctx['client'] = client
    ctx['token_owner_w1'] = fetch_token(client, ctx['owner_w1_email'], ctx['owner_w1_pw'])
    ctx['token_owner_w2'] = fetch_token(client, ctx['owner_w2_email'], ctx['owner_w2_pw'])
    ctx['token_viewer_w1'] = fetch_token(client, ctx['viewer_w1_email'], ctx['viewer_w1_pw'])
    ctx['token_owner_both'] = fetch_token(client, ctx['owner_both_email'], ctx['owner_both_pw'])
    return ctx
