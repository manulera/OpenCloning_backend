"""Smoke tests for the Typer command surface."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from opencloning_cli.main import app


runner = CliRunner()


def _invoke(*args: str):
    # mix_stderr=False would be nicer on Click>=8.2 but Typer's CliRunner
    # surfaces both streams through result.output, which is enough for us.
    return runner.invoke(app, list(args))


class TestHelpAndTree:
    def test_root_help(self):
        result = _invoke('--help')
        assert result.exit_code == 0
        assert 'db' in result.output

    def test_nested_help(self):
        result = _invoke('db', 'test', 'snapshot', '--help')
        assert result.exit_code == 0
        assert 'create' in result.output
        assert 'restore' in result.output


class TestSeedCommand:
    def test_success_human_output(self, temp_workspace):
        _, config = temp_workspace

        result = _invoke('db', 'test', 'seed')

        assert result.exit_code == 0, result.output
        assert result.output.strip() == ''
        assert Path(config.database_path).exists()


class TestSnapshotCommands:
    def test_create_and_restore(self, temp_workspace):
        _, _ = temp_workspace

        assert _invoke('db', 'test', 'seed').exit_code == 0
        create = _invoke('db', 'test', 'snapshot', 'create')
        assert create.exit_code == 0, create.output
        snapshot_dir = Path(temp_workspace[1].database_path).parent / 'snapshot'
        assert (snapshot_dir / 'db' / Path(temp_workspace[1].database_path).name).exists()

        restore = _invoke('db', 'test', 'snapshot', 'restore')
        assert restore.exit_code == 0, restore.output

    def test_restore_missing_snapshot_exit_code(self, temp_workspace):
        result = _invoke('db', 'test', 'snapshot', 'restore')

        assert result.exit_code != 0

    def test_custom_snapshot_dir(self, temp_workspace, tmp_path):
        _, _ = temp_workspace
        assert _invoke('db', 'test', 'seed').exit_code == 0

        custom = tmp_path / 'custom-snapshot'
        result = _invoke('db', 'test', 'snapshot', 'create', '--snapshot-dir', str(custom))

        assert result.exit_code == 0
        assert (custom / 'db' / Path(temp_workspace[1].database_path).name).exists()


class TestResetCommand:
    def test_first_invocation_seeds(self, temp_workspace):
        _, config = temp_workspace

        result = _invoke('db', 'test', 'reset')

        assert result.exit_code == 0, result.output
        assert result.output.strip() == ''
        assert Path(config.database_path).exists()

    def test_second_invocation_restores(self, temp_workspace):
        first = _invoke('db', 'test', 'reset')
        assert first.exit_code == 0

        second = _invoke('db', 'test', 'reset')
        assert second.exit_code == 0
        assert second.output.strip() == ''
