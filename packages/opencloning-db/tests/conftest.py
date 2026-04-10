from typing import Generator
from sqlalchemy.orm import Session

from opencloning_db.config import Config, get_config, set_config
import opencloning_db.db as db_module
from fastapi.testclient import TestClient
from opencloning_db.api import app
from opencloning_db.deps import get_db
from opencloning_db.models import Base
from sqlalchemy.engine import Engine
import tempfile
import pytest

_JWT_SECRET = 'test-jwt-secret-not-for-production'


@pytest.fixture
def engine_client_config() -> Generator[tuple[Engine, TestClient, Config], None, None]:
    """Temp SQLite DB, upload dirs, ``get_db`` override, and ``TestClient``.

    Also use via ``@pytest.mark.usefixtures("engine_client_config")`` when tests
    only need ``get_config()`` paths (e.g. model tests with their own engine).
    """
    default_config = get_config()
    with tempfile.TemporaryDirectory() as tmp_dir:
        test_config = Config(
            database_url=f'sqlite:///{tmp_dir}/test.db',
            jwt_secret=_JWT_SECRET,
            sequence_files_dir=tmp_dir,
            sequencing_files_dir=tmp_dir,
        )
        set_config(test_config)
        engine = db_module.get_engine(test_config)
        Base.metadata.create_all(engine)

        def override_get_db():
            session = Session(engine)
            try:
                yield session
            finally:
                session.close()

        app.dependency_overrides[get_db] = override_get_db
        yield engine, TestClient(app), test_config
        app.dependency_overrides.pop(get_db, None)
        set_config(default_config)
