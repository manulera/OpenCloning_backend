"""
Initialize the database: create tables and seed from cloning_strategy.json.
Run this script manually when you need to (re)create the database.
"""

import json
import os
import glob
import opencloning_linkml.datamodel.models as opencloning_models
from sqlalchemy.orm import Session
from sqlalchemy.schema import CreateTable
from sqlalchemy.dialects import postgresql

from auth.security import get_password_hash
from config import Config, get_config
from models import (
    Base,
    Line,
    Primer,
    SequenceInLine,
    SequenceSample,
    Tag,
    SequenceType,
    User,
    Workspace,
    WorkspaceMembership,
    WorkspaceRole,
)
from db import cloning_strategy_to_db, create_sequencing_file, get_engine


def init_db(config: Config):
    if config.database_path is not None and os.path.exists(config.database_path):
        os.remove(config.database_path)

    os.makedirs(config.sequence_files_dir, exist_ok=True)
    os.makedirs(config.sequencing_files_dir, exist_ok=True)

    engine = get_engine(config)
    Base.metadata.create_all(engine)

    cloning_strategies = []
    file_names = []
    for file in glob.glob('init_db/*.json'):
        with open(file) as f:
            cloning_strategies.append(opencloning_models.CloningStrategy.model_validate(json.load(f)))
            file_names.append(os.path.basename(file))

    with Session(engine) as session:
        last_seq = None
        sequencing_sequence = None
        # Dev-only: replace with env-driven or unset password before production.
        bootstrap_user = User(
            email='bootstrap@example.com',
            display_name='Bootstrap User',
            password_hash=get_password_hash('password'),
            is_instance_admin=True,
        )
        workspace = Workspace(name='Bootstrap Workspace')
        session.add_all([bootstrap_user, workspace])
        session.flush()
        session.add(
            WorkspaceMembership(
                user_id=bootstrap_user.id,
                workspace_id=workspace.id,
                role=WorkspaceRole.owner,
            )
        )

        parent_strain = Line(uid='parent_strain', workspace_id=workspace.id)
        session.add(parent_strain)
        for cloning_strategy, file_name in zip(cloning_strategies, file_names):
            tag_name = os.path.basename(file_name).split('.')[0]
            tag = Tag(name=tag_name, workspace_id=workspace.id)
            session.add(tag)
            sequences, id_mappings = cloning_strategy_to_db(cloning_strategy, session, workspace.id)
            new_line = Line(uid=f"{tag_name}-line", workspace_id=workspace.id, parents=[parent_strain])
            for seq in sequences:
                if seq.sequence_type == SequenceType.allele or seq.sequence_type == SequenceType.plasmid:
                    new_line.sequences_in_line.append(SequenceInLine(sequence=seq))
            new_line.tags.append(tag)
            session.add(new_line)

            for seq in sequences:
                seq.tags.append(tag)
            last_seq = sequences[-1]
            if cloning_strategy.description == 'sequencing':
                sequencing_sequence = last_seq
            session.add(
                SequenceSample(
                    uid=f"{tag_name}-sample",
                    sequence_id=last_seq.id,
                    uid_workspace_id=workspace.id,
                )
            )

        # Take the last 4 primers and give them a uid
        for primer in session.query(Primer).order_by(Primer.id).limit(4).all():
            primer.uid = f"ML{primer.id}"
            primer.uid_workspace_id = workspace.id
            primer.tags.append(tag)
            session.add(primer)

        for file in glob.glob('init_db/sequencing_data/*'):
            with open(file, 'rb') as f:
                content = f.read()
            file_name = os.path.basename(file)
            session.add(create_sequencing_file(sequencing_sequence, content, file_name))
        session.commit()

    print('Database initialized successfully.')

    sql_file_content = ''

    # Print DDL for all tables in the database
    for table in Base.metadata.tables.values():
        sql_file_content += str(CreateTable(table).compile(dialect=postgresql.dialect()))

    with open('postgres_ddl.sql', 'w') as f:
        f.write(sql_file_content)


if __name__ == '__main__':  # pragma: no cover
    config = get_config()
    init_db(config)
