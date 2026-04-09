from init_db import init_db
from config import Config, set_config
import tempfile
import os


def test_init_db():
    with tempfile.TemporaryDirectory() as tmp_dir_sequences:
        with tempfile.TemporaryDirectory() as tmp_dir_sequencing:
            with tempfile.NamedTemporaryFile() as tmp_db_file:
                config = Config(
                    database_url=f"sqlite:///{tmp_db_file.name}",
                    sequence_files_dir=tmp_dir_sequences,
                    sequencing_files_dir=tmp_dir_sequencing,
                )
                set_config(config)

                init_db(config)
                assert len(os.listdir(tmp_dir_sequences)) == 48
                assert len(os.listdir(tmp_dir_sequencing)) == 3
