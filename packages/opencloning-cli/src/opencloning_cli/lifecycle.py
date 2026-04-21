"""DB test lifecycle primitives.

These functions are pure library code (no Typer or click imports) so they can
be unit-tested directly against a ``Config`` pointing at a temporary
workspace.
"""

from __future__ import annotations

import io
import shutil
from contextlib import redirect_stdout
from pathlib import Path

import opencloning_db.db as _db_module
from opencloning_db.config import Config
from opencloning_db.init_db import init_db as _init_db

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
