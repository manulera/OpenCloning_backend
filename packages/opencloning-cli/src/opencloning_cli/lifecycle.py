"""DB test lifecycle primitives.

These functions are pure library code (no Typer or click imports) so they can
be unit-tested directly against a ``Config`` pointing at a temporary
workspace.
"""

from __future__ import annotations

import io
import json
import base64
import shutil
from contextlib import redirect_stdout
from pathlib import Path
from typing import Any

import opencloning_db.db as _db_module
from opencloning_db.config import Config, get_config
from opencloning_db.init_db import init_db as _init_db
from .stubs import stubs

# Subdirectory names inside a snapshot directory.
_DB_SUBDIR = 'db'
_SEQUENCE_SUBDIR = 'sequence_files'
_SEQUENCING_SUBDIR = 'sequencing_files'

_DEFAULT_SNAPSHOT_SUBDIR = 'snapshot'


class SnapshotMissingError(Exception):
    """Raised when snapshot files are missing or incomplete."""


def default_snapshot_dir(config: Config) -> Path:
    """Return the default snapshot directory for *config*.

    Derived from the parent of the SQLite database file to keep all test data
    colocated under a single ``OPENCLONING_TEST_DATA_ROOT``.
    """
    db_path = config.database_path
    if db_path is None:
        raise ValueError(
            'opencloning-cli db test commands require a file-based SQLite '
            'database_url (sqlite:///...). Non-file backends are not supported.'
        )
    return Path(db_path).expanduser().parent / _DEFAULT_SNAPSHOT_SUBDIR


def resolve_snapshot_dir(config: Config, override: Path | None) -> Path:
    """Return *override* when provided, otherwise the default snapshot dir."""
    if override is not None:
        return Path(override).expanduser()
    return default_snapshot_dir(config)


def _dispose_engine() -> None:
    """Dispose any cached SQLAlchemy engine so DB files can be replaced.

    ``opencloning_db.db`` caches a module-level engine keyed by URL; leaving
    it open would keep file handles against the live SQLite file on Windows
    and can cause WAL sidecar files to linger on POSIX.
    """
    if _db_module._engine is not None:
        _db_module._engine.dispose()
        _db_module._engine = None
        _db_module._bound_database_url = None


def _sqlite_files(db_path: Path) -> list[Path]:
    """Return the SQLite file and any live WAL/SHM sidecars for *db_path*."""
    candidates = [db_path, db_path.with_name(db_path.name + '-wal'), db_path.with_name(db_path.name + '-shm')]
    return [p for p in candidates if p.exists()]


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _copy_tree(src: Path, dst: Path) -> None:
    """Copy *src* directory onto *dst*, replacing *dst* if it exists."""
    if dst.exists():
        shutil.rmtree(dst)
    if src.exists():
        shutil.copytree(src, dst)
    else:
        dst.mkdir(parents=True, exist_ok=True)


def seed(config: Config) -> None:
    """Run ``opencloning_db.init_db.init_db`` against *config*.

    Removes any existing SQLite DB (``init_db`` does this internally) and
    writes a fresh seeded database plus the sequence/sequencing file dirs.
    """
    _dispose_engine()
    db_path_str = config.database_path
    if db_path_str is not None:
        Path(db_path_str).expanduser().parent.mkdir(parents=True, exist_ok=True)
    # ``init_db`` prints a success message; keep CLI successful runs silent.
    with redirect_stdout(io.StringIO()):
        _init_db(config)
    # Dispose again so the next caller sees a fresh engine bound to the
    # newly-created DB file rather than a stale handle from init_db.
    _dispose_engine()


def snapshot_create(config: Config, snapshot_dir: Path) -> None:
    """Capture the current DB and file directories into *snapshot_dir*.

    The resulting directory layout is::

        <snapshot_dir>/
          db/<db_basename>(+ optional -wal/-shm)
          sequence_files/...
          sequencing_files/...
    """
    db_path_str = config.database_path
    if db_path_str is None:
        raise ValueError('snapshot_create requires a file-based SQLite database_url (sqlite:///...).')
    db_path = Path(db_path_str)
    if not db_path.exists():
        raise FileNotFoundError(
            f'Cannot snapshot: database file does not exist at {db_path}. ' 'Run "opencloning-cli db test seed" first.'
        )

    _dispose_engine()

    snapshot_dir = Path(snapshot_dir)
    if snapshot_dir.exists():
        shutil.rmtree(snapshot_dir)
    snapshot_dir.mkdir(parents=True)

    db_target_dir = snapshot_dir / _DB_SUBDIR
    db_target_dir.mkdir()
    for live_file in _sqlite_files(db_path):
        shutil.copy2(live_file, db_target_dir / live_file.name)

    _copy_tree(Path(config.sequence_files_dir), snapshot_dir / _SEQUENCE_SUBDIR)
    _copy_tree(Path(config.sequencing_files_dir), snapshot_dir / _SEQUENCING_SUBDIR)


def snapshot_restore(config: Config, snapshot_dir: Path) -> None:
    """Replace live DB and file dirs with the contents of *snapshot_dir*.

    Verifies snapshot files before touching live data. Disposes the cached
    SQLAlchemy engine so callers get a fresh handle against the restored DB.
    """
    snapshot_dir = Path(snapshot_dir)
    if not snapshot_dir.exists():
        raise SnapshotMissingError(f'Snapshot directory does not exist: {snapshot_dir}.')

    db_path_str = config.database_path
    if db_path_str is None:
        raise ValueError('snapshot_restore requires a file-based SQLite database_url (sqlite:///...).')
    db_path = Path(db_path_str)

    snapshot_db_dir = snapshot_dir / _DB_SUBDIR
    snapshot_db_file = snapshot_db_dir / db_path.name
    if not snapshot_db_file.exists():
        raise SnapshotMissingError(f'Snapshot is missing its DB file at {snapshot_db_file}.')

    _dispose_engine()

    for live_file in _sqlite_files(db_path):
        live_file.unlink()

    _ensure_parent(db_path)
    for src in snapshot_db_dir.iterdir():
        shutil.copy2(src, db_path.parent / src.name)

    _copy_tree(snapshot_dir / _SEQUENCE_SUBDIR, Path(config.sequence_files_dir))
    _copy_tree(snapshot_dir / _SEQUENCING_SUBDIR, Path(config.sequencing_files_dir))


def reset(config: Config, snapshot_dir: Path) -> None:
    """Fast-path: restore from *snapshot_dir*; fall back to reseed+snapshot.

    When the snapshot is missing or its manifest is corrupt, this function
    runs :func:`seed` and then :func:`snapshot_create` so subsequent resets
    take the fast path.
    """
    snapshot_dir = Path(snapshot_dir)
    try:
        snapshot_restore(config, snapshot_dir)
    except SnapshotMissingError:
        seed(config)
        snapshot_create(config, snapshot_dir)


def _sanitize_headers(headers: dict[str, str] | None) -> dict[str, str]:
    """Drop noisy values and replace auth tokens with placeholders."""
    if not headers:
        return {}

    sanitized: dict[str, str] = {}
    for key, value in headers.items():
        normalized = key.lower()
        if normalized in {'authorization'}:
            sanitized[normalized] = 'Bearer __TEST_TOKEN__'
            continue
        if normalized in {'x-workspace-id', 'content-type'}:
            sanitized[normalized] = value
    return sanitized


def create_stub(
    test_client: Any,
    endpoint: str,
    method: str,
    params: dict[str, Any] | None = None,
    body: dict[str, Any] | list[Any] | None = None,
    headers: dict[str, str] | None = None,
    multipart_files: list[dict[str, str]] | None = None,
    binary_response: bool = False,
) -> dict[str, Any]:
    """Perform one request through *test_client* and return a stub payload."""

    method_name = method.lower()
    requester = getattr(test_client, method_name)
    request_headers = headers or {}

    request_kwargs: dict[str, Any] = {'params': params, 'headers': request_headers}
    if body is not None:
        request_kwargs['json'] = body
    if multipart_files:
        files: list[tuple[str, tuple[str, bytes, str]]] = []
        for file_spec in multipart_files:
            files.append(
                (
                    'files',
                    (
                        file_spec['filename'],
                        file_spec['content'].encode('utf-8'),
                        file_spec.get('content_type', 'application/octet-stream'),
                    ),
                )
            )
        request_kwargs['files'] = files

    response = requester(endpoint, **request_kwargs)

    if binary_response:
        response_body = base64.b64encode(response.content).decode('ascii')
    else:
        try:
            response_body = response.json()
        except ValueError:
            response_body = response.text

    request_body: Any = body
    if multipart_files:
        request_body = {
            'multipart_files': [
                {
                    'filename': file_spec['filename'],
                    'content_type': file_spec.get('content_type', 'application/octet-stream'),
                    'content': file_spec['content'],
                }
                for file_spec in multipart_files
            ]
        }

    payload = {
        'method': method.upper(),
        'endpoint': endpoint,
        'params': params,
        'headers': _sanitize_headers(request_headers),
        'body': request_body,
        'response': {
            'status_code': response.status_code,
            'headers': _sanitize_headers(dict(response.headers)),
            'body': response_body,
        },
    }
    if binary_response:
        payload['response']['body_encoding'] = 'base64'
    return payload


def _default_auth_headers(test_client: Any) -> dict[str, str]:
    token_response = test_client.post(
        '/auth/token',
        data={'username': 'bootstrap@example.com', 'password': 'password'},
    )
    token_response.raise_for_status()
    token = token_response.json()['access_token']

    workspaces_response = test_client.get(
        '/workspaces',
        headers={'Authorization': f'Bearer {token}'},
    )
    workspaces_response.raise_for_status()
    workspace_id = workspaces_response.json()[0]['id']

    return {
        'Authorization': f'Bearer {token}',
        'X-Workspace-Id': str(workspace_id),
    }


def _example_body(name: str) -> dict[str, Any]:
    if name == 'cs_pcr':
        from Bio.Seq import reverse_complement
        from pydna import opencloning_models
        from pydna.assembly2 import pcr_assembly
        from pydna.dseqrecord import Dseqrecord
        from pydna.primer import Primer

        primer1 = Primer('ACGTACGT')
        primer2 = Primer(reverse_complement('GCGCGCGC'))
        pcr_template = Dseqrecord('ccccACGTACGTAAAAAAGCGCGCGCcccc', circular=True)
        pcr_product, *_ = pcr_assembly(pcr_template, primer1, primer2, limit=8)
        cs_pcr = opencloning_models.CloningStrategy.from_dseqrecords([pcr_product])
        return opencloning_models.CloningStrategy.model_validate(cs_pcr.model_dump(mode='json')).model_dump(
            mode='json'
        )
    raise ValueError(f'Unknown example body "{name}".')


def write_stubs(output_dir: Path):
    """Generate and persist one predefined DB test stub JSON."""
    try:
        from fastapi.testclient import TestClient
        from opencloning_db.api import app
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            'Generating stubs requires opencloning-db and fastapi test dependencies installed.'
        ) from exc

    config = get_config()
    reset(config, resolve_snapshot_dir(config, None))

    target_dir = Path(output_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    client = TestClient(app)
    headers = _default_auth_headers(client)
    generated_payloads: dict[str, dict[str, Any]] = {}
    runtime_vars: dict[str, Any] = {}
    for stub in stubs:
        output_file = target_dir / f'{stub.name}.json'
        try:
            endpoint = stub.endpoint.format(**runtime_vars)
            body = stub.body
            if stub.body_from_stub:
                source_payload = generated_payloads.get(stub.body_from_stub)
                if source_payload is None:
                    raise ValueError(f'Unknown body source stub "{stub.body_from_stub}" for "{stub.name}".')
                body = source_payload['response']['body']
            if stub.body_from_example:
                body = _example_body(stub.body_from_example)

            payload = create_stub(
                client,
                endpoint=endpoint,
                method=stub.method,
                headers=headers,
                params=stub.params,
                body=body,
                multipart_files=stub.multipart_files,
                binary_response=stub.binary_response,
            )
            generated_payloads[stub.name] = payload

            response_body = payload.get('response', {}).get('body')
            if isinstance(response_body, dict) and 'id' in response_body:
                runtime_vars['last_id'] = response_body['id']
            if isinstance(response_body, list) and response_body and isinstance(response_body[0], dict):
                first_id = response_body[0].get('id')
                if first_id is not None:
                    runtime_vars['last_file_id'] = first_id

            with output_file.open('w', encoding='utf-8') as handle:
                json.dump(payload, handle, indent=2, sort_keys=True)
                handle.write('\n')
        finally:
            if stub.reset_db:
                reset(config, resolve_snapshot_dir(config, None))
