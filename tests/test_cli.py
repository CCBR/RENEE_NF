import subprocess
from unittest.mock import Mock

from click.testing import CliRunner

from src.__main__ import cli, repo_base


def invoke(*args):
    """Invoke the CLI using the installed command's program name."""
    result = CliRunner().invoke(cli, list(args), prog_name="renee")
    assert result.exit_code == 0, result.output
    return result.output


def test_help():
    output = invoke("--help")
    assert "RENEE" in output
    assert "Usage: renee [OPTIONS]" in output
    assert "renee_nf" not in output


def test_version():
    output = invoke("--version")
    assert f"renee, version {repo_base('VERSION').read_text().strip()}" in output


def test_installed_command():
    """Verify the installed console script uses the renamed package entry point."""
    result = subprocess.run(
        ["renee", "--version"], capture_output=True, text=True, check=True
    )
    assert f"renee, version {repo_base('VERSION').read_text().strip()}" in result.stdout


def test_citation():
    output = invoke("--citation")
    assert "@misc{" in output
    assert "RENEE:" in output
    assert "RENEE_NF:" not in output


def test_subcommands_help():
    for cmd in ["run", "init", "cite"]:
        assert f"renee {cmd} [OPTIONS]" in invoke(cmd, "--help")


def test_run_repository_name(tmp_path, monkeypatch):
    """Accept the current repository URL and forward the RENEE pipeline name."""
    (tmp_path / "nextflow.config").touch()
    nextflow_run = Mock()
    monkeypatch.setattr("ccbr_tools.pipeline.nextflow.run", nextflow_run)

    invoke(
        "run", "--main", "CCBR/RENEE_NF", "--output", str(tmp_path), "--mode", "local"
    )

    nextflow_run.assert_called_once_with(
        nextfile_path="CCBR/RENEE_NF",
        mode="local",
        force_all=False,
        pipeline_name="RENEE",
        hpc_memory="4GB",
        nextflow_args=(),
    )
