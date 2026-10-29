import csv
import json
import os
import re
import shlex
import subprocess
from pathlib import Path
from unittest.mock import Mock

import click
import pytest
from click.testing import CliRunner

from src.__main__ import cli, repo_base
from src.launcher import cache_images, execute, prepare_args, write_samplesheet


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
    for cmd in ["run", "init", "cite", "build", "cache", "debug", "unlock"]:
        assert f"renee {cmd} [OPTIONS]" in invoke(cmd, "--help")


def test_run_repository_name(tmp_path, monkeypatch):
    """Accept the current repository URL and forward the RENEE pipeline name."""
    (tmp_path / "nextflow.config").touch()
    nextflow_run = Mock()
    monkeypatch.setattr("src.__main__.execute", nextflow_run)

    invoke(
        "run", "--main", "CCBR/RENEE_NF", "--output", str(tmp_path), "--mode", "local"
    )

    nextflow_run.assert_called_once_with(
        "CCBR/RENEE_NF", tmp_path, "local", False, [], False, None, None, None
    )


def test_auto_init_and_dry_run(tmp_path, monkeypatch):
    """A first run initializes the directory and previews without submission."""
    launch = Mock()
    monkeypatch.setattr("src.__main__.execute", launch)
    output = tmp_path / "new run"
    invoke("run", "--output", str(output), "--dry-run", "--mode", "local")
    assert (output / "nextflow.config").is_file()
    assert launch.call_args.args[4] == ["-preview"]
    config = output / "nextflow.config"
    config.write_text("// user configuration\n")
    invoke("run", "--output", str(output), "--runmode", "init")
    assert config.read_text() == "// user configuration\n"
    launch.assert_called_once()


def test_build_forwarding(tmp_path, monkeypatch):
    """Build preserves native reference parameters and driver options."""
    launch = Mock()
    monkeypatch.setattr("src.__main__.execute", launch)
    invoke(
        "build",
        "--output",
        str(tmp_path),
        "--mode",
        "local",
        "--ref_fa",
        str(repo_base("tests/test_data/hg38_stub/genome.fa")),
        "--ref_gtf",
        str(repo_base("tests/test_data/hg38_stub/genes.gtf")),
        "--ref_name",
        "hg38",
        "--gtf_ver",
        "47",
        "--small_genome",
        "--partition",
        "student",
        "--time",
        "12:00:00",
        "--wait",
    )
    args = launch.call_args.args
    assert args[4][0] == "--build"
    assert "--small_genome" in args[4]
    assert args[5:8] == (True, "student", "12:00:00")


@pytest.mark.parametrize("paired", [True, False])
def test_fastq_inputs(tmp_path, monkeypatch, paired):
    """Legacy FASTQ lists are converted without losing the sample name."""
    files = [tmp_path / "sample.R1.fastq.gz"]
    if paired:
        files.append(tmp_path / "sample.R2.fastq.gz")
    for fastq in files:
        fastq.touch()
    output = tmp_path / "output"
    (output / "log").mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    args = prepare_args(
        ("--input", "sample.R?.fastq.gz", "--small_rna", "-profile", "ci_stub"),
        output,
        repo_base(),
    )
    with Path(args[1]).open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows == [
        {
            "sample": "sample",
            "replicate": "",
            "fastq_1": str(files[0]),
            "fastq_2": str(files[1]) if paired else "",
        }
    ]
    assert args[2:] == ["--small_rna", "-profile", "ci_stub"]


@pytest.mark.parametrize(
    "names, message",
    [
        (["a.R2.fastq.gz"], "Missing FASTQ mate"),
        (["a.R1.fastq.gz", "b.R1.fastq.gz", "b.R2.fastq.gz"], "Missing FASTQ mate"),
        (["a.R1.fastq.gz", "a.R1.fastq.gz"], "Duplicate FASTQ"),
        (["a.fastq.gz"], "Expected an existing"),
    ],
)
def test_invalid_fastqs(tmp_path, names, message):
    """Fail before launch rather than silently dropping invalid samples."""
    files = [tmp_path / name for name in names]
    for fastq in files:
        fastq.touch()
    with pytest.raises(click.BadParameter, match=message):
        write_samplesheet(files, tmp_path / "samples.csv")


def test_cli_paths_and_native_arguments(tmp_path, monkeypatch):
    """Resolve paths at invocation, but preserve quoted values and repeated configs."""
    monkeypatch.chdir(tmp_path)
    args = prepare_args(
        (
            "--ref_fa=genome.fa",
            "-c",
            "one.config",
            "-c",
            "two.config",
            "-params-file",
            "params.yml",
            "--star_pass2_args",
            "--foo bar",
        ),
        tmp_path,
        repo_base(),
    )
    assert args == [
        "--ref_fa",
        str(tmp_path / "genome.fa"),
        "-c",
        str(tmp_path / "one.config"),
        "-c",
        str(tmp_path / "two.config"),
        "-params-file",
        str(tmp_path / "params.yml"),
        "--star_pass2_args",
        "--foo bar",
    ]


def test_legacy_genome_json(tmp_path):
    """Import the legacy schema without rebuilding or mutating the source file."""
    source = tmp_path / "custom.json"
    source.write_text(
        json.dumps(
            {
                "references": {
                    "rnaseq": {
                        "GENOMEFILE": "/references/genome.fa",
                        "GTFFILE": "/references/genes.gtf",
                        "GENEINFO": "/references/geneinfo.bed",
                        "RSEMREF": "/references/rsem/prefix",
                        "ORGANISM": "custom",
                    }
                }
            }
        )
    )
    (tmp_path / "log").mkdir()
    args = prepare_args(("--genome", str(source)), tmp_path, repo_base())
    assert args[0] == "-c"
    assert args[2:] == ["--genome", "custom"]
    config = Path(args[1]).read_text()
    assert "params.genomes['custom']" in config
    assert "rsem_ref: '/references/rsem/prefix'" in config


def test_local_quoting_and_failures(tmp_path, monkeypatch):
    """Execute a fake Nextflow binary to verify real shell argument boundaries."""
    executable = tmp_path / "nextflow"
    executable.write_text(
        "#!/usr/bin/env python\nimport json, sys\n"
        "from pathlib import Path\n"
        "Path('argv.json').write_text(json.dumps(sys.argv[1:]))\n"
        "sys.exit(9 if '--fail' in sys.argv else 0)\n"
    )
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setattr("src.launcher.get_hpc", lambda: None)
    execute(
        "main file.nf",
        tmp_path,
        "local",
        False,
        ["--star_pass2_args", "--foo bar", "-c", "one config", "-c", "two config"],
    )
    args = json.loads((tmp_path / "argv.json").read_text())
    assert args == [
        "run",
        "main file.nf",
        "--star_pass2_args",
        "--foo bar",
        "-c",
        "one config",
        "-c",
        "two config",
        "-resume",
    ]
    with pytest.raises(click.ClickException, match="Workflow launch failed"):
        execute("main.nf", tmp_path, "local", False, ["--fail"])


def test_slurm_wait_and_preview(tmp_path, monkeypatch):
    """Wait, partition and walltime reach sbatch; previews stay off the queue."""
    monkeypatch.setattr("src.launcher.get_hpc", lambda: None)
    process = Mock()
    monkeypatch.setattr("src.launcher.subprocess.run", process)
    execute(
        "main.nf",
        tmp_path,
        "slurm",
        True,
        ["-profile", "singularity"],
        wait=True,
        partition="student",
        walltime="12:00:00",
    )
    command = process.call_args.args[0]
    assert "--wait" in command
    assert "--partition=student" in command
    assert "--time=12:00:00" in command
    script = (tmp_path / "submit_slurm.sh").read_text()
    assert "-resume" not in script
    assert "singularity,slurm" in script
    process.reset_mock()
    execute("main.nf", tmp_path, "slurm", False, ["-preview"])
    assert process.call_args.args[0][:2] == ["bash", "-c"]


def test_private_cache_and_resume(tmp_path, monkeypatch):
    """Private cache overrides propagate to the process without changing parent env."""
    monkeypatch.setattr("src.launcher.get_hpc", lambda: None)
    process = Mock()
    monkeypatch.setattr("src.launcher.subprocess.run", process)
    execute(
        "main.nf",
        tmp_path,
        "local",
        True,
        ["-resume=old"],
        singularity_cache=str(tmp_path / "private cache"),
    )
    script = process.call_args.args[0][2]
    assert (
        f"export SINGULARITY_CACHEDIR={shlex.quote(str(tmp_path / 'private cache'))}"
        in script
    )
    assert "-resume" not in script


def test_cache_inspection(tmp_path, monkeypatch):
    """Use Nextflow cache names and avoid pulls or directories for dry-run."""
    process = Mock(
        return_value=Mock(
            stdout=json.dumps(
                {
                    "processes": [
                        {"container": "quay.io/org/tool:1"},
                        {"container": "quay.io/org/tool:1"},
                    ]
                }
            )
        )
    )
    monkeypatch.setattr("src.launcher.subprocess.run", process)
    destination = tmp_path / "cache"
    cache_images(repo_base(), destination, True)
    assert not destination.exists()
    process.assert_called_once()
    cache_images(repo_base(), destination, False)
    assert process.call_args.args[0] == [
        "singularity",
        "pull",
        str(destination / "quay.io-org-tool-1.img"),
        "docker://quay.io/org/tool:1",
    ]


def test_invalid_mode(tmp_path):
    """Reject invalid modes before initializing or running a workflow."""
    result = CliRunner().invoke(
        cli, ["run", "--mode", "invalid", "--output", str(tmp_path)]
    )
    assert result.exit_code == 2
    assert "Invalid value" in result.output


def test_generated_samplesheet_schema_accepts_spaces(tmp_path):
    """Generated CSV paths remain valid when the launch directory contains spaces."""
    schema = json.loads(repo_base("nextflow_schema.json").read_text())
    pattern = schema["$defs"]["input_output_options"]["properties"]["input"]["pattern"]
    assert re.fullmatch(
        pattern, str(tmp_path / "my analysis/log/launcher_samplesheet.csv")
    )


def test_schema_definitions_are_referenced_at_root():
    """Keep every parameter group in the root allOf required by nf-core lint."""
    schema = json.loads(repo_base("nextflow_schema.json").read_text())
    references = [entry["$ref"] for entry in schema["allOf"]]
    assert sorted(references) == sorted(f"#/$defs/{name}" for name in schema["$defs"])
