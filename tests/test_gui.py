"""Tests for the original PySimpleGUI form's Nextflow execution adapter."""

import subprocess
import sys
from unittest.mock import Mock, call

import click
import pytest
from click.testing import CliRunner

from src.__main__ import cli, repo_base
from src.gui import genome_choices, get_fastqs, launch_gui, run_nextflow


@pytest.fixture
def sg(monkeypatch):
    """Provide deterministic PySimpleGUI events without requiring a display."""
    module = Mock()
    module.WIN_CLOSED = "__CLOSED__"
    monkeypatch.setitem(sys.modules, "PySimpleGUI", module)
    return module


@pytest.fixture
def selections(tmp_path):
    """Use legacy FASTQ names and output paths containing spaces."""
    reads = tmp_path / "input reads"
    reads.mkdir()
    (reads / "sample.R1.fastq.gz").touch()
    return {
        "--INDIR--": str(reads),
        "--OUTDIR--": str(tmp_path / "output folder"),
        "--ANNOTATION--": "hg38_36",
    }


def cancel_after(sg, events):
    """Finish an event sequence with the original Cancel action."""
    sg.Window.return_value.read.side_effect = [*events, ("--CANCEL--", {})]
    with pytest.raises(SystemExit) as error:
        launch_gui(repo_base(), debug=False)
    assert error.value.code == 69


def test_original_form(sg):
    """Keep the original fields, buttons, window location, and close status."""
    sg.Window.return_value.read.return_value = (sg.WIN_CLOSED, {})
    with pytest.raises(SystemExit, match="69"):
        launch_gui(repo_base(), debug=False)
    assert sg.InputText.call_args_list == [
        call(key="--INDIR--"),
        call(key="--OUTDIR--"),
    ]
    assert sg.FolderBrowse.call_args_list == [
        call(target="--INDIR--"),
        call(target="--OUTDIR--"),
    ]
    assert sg.Combo.call_args.kwargs["key"] == "--ANNOTATION--"
    assert sg.Window.call_args.kwargs == {"location": (0, 500), "finalize": True}
    sg.Submit.assert_called_once_with(key="--SUBMIT--", font=("Helvetica", 12))
    sg.Cancel.assert_called_once_with(key="--CANCEL--", font=("Helvetica", 12))


def test_preview_then_confirmed_submission(sg, selections, monkeypatch):
    """Preview and submit the same selections through the Nextflow CLI."""
    execute = Mock(
        return_value=subprocess.CompletedProcess([], 0, "HANDLE_PRESEQ_ERROR")
    )
    monkeypatch.setattr("src.gui.run_nextflow", execute)
    sg.Window.return_value.read.return_value = ("--SUBMIT--", selections)
    sg.popup_yes_no.return_value = "Yes"
    launch_gui(repo_base(), debug=False)
    preview = execute.call_args_list[0].args[0]
    assert preview == [
        sys.executable,
        str(repo_base("main.py")),
        "run",
        "--input",
        *get_fastqs(selections["--INDIR--"]),
        "--output",
        selections["--OUTDIR--"],
        "--genome",
        "hg38_36",
        "--mode",
        "slurm",
        "--threads",
        "2",
        "--dry-run",
    ]
    assert execute.call_args_list[1] == call(preview[:-1])
    assert execute.call_count == 2
    assert sg.popup_scrolled.call_args_list[0].kwargs["title"] == "Dryrun:STDOUT/STDERR"
    sg.popup_yes_no.assert_called_once()
    assert sg.popup_yes_no.call_args.args == ("Submit run to slurm?",)
    sg.Window.return_value.close.assert_called_once()


def test_declined_submission_clears_form(sg, selections, monkeypatch):
    """Retain the original No action rather than adding a second submit button."""
    execute = Mock(return_value=subprocess.CompletedProcess([], 0, "Preview"))
    monkeypatch.setattr("src.gui.run_nextflow", execute)
    sg.popup_yes_no.return_value = "No"
    cancel_after(sg, [("--SUBMIT--", selections)])
    execute.assert_called_once()
    assert sg.Window.return_value.Element.call_args_list == [
        call("--INDIR--"),
        call("--OUTDIR--"),
        call("--ANNOTATION--"),
    ]
    assert sg.Window.return_value.Element.return_value.update.call_args_list == [
        call(""),
        call(""),
        call(""),
    ]


@pytest.mark.parametrize("output", ["", "ERROR: invalid reference"])
def test_failed_preview_never_submits(sg, selections, monkeypatch, output):
    """Honor launcher failures even when they have no output."""
    execute = Mock(return_value=subprocess.CompletedProcess([], 1, output))
    monkeypatch.setattr("src.gui.run_nextflow", execute)
    cancel_after(sg, [("--SUBMIT--", selections)])
    execute.assert_called_once()
    sg.popup_yes_no.assert_not_called()


def test_submission_failure(sg, selections, monkeypatch):
    """Do not show a successful goodbye when Nextflow submission fails."""
    execute = Mock(
        side_effect=[
            subprocess.CompletedProcess([], 0, "Preview"),
            subprocess.CompletedProcess([], 1, "sbatch failed"),
        ]
    )
    monkeypatch.setattr("src.gui.run_nextflow", execute)
    sg.popup_yes_no.return_value = "Yes"
    cancel_after(sg, [("--SUBMIT--", selections)])
    assert sg.PopupError.call_args.args == (
        "Workflow submission failed (exit status 1).",
    )
    assert execute.call_count == 2
    # Only the later Cancel action displays the goodbye.
    sg.popup_auto_close.assert_called_once()


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("--INDIR--", "", "Input folder must be provided!!"),
        ("--INDIR--", "/nonexistent/renee-test-reads", "Input folder doesn't exist!!"),
        ("--OUTDIR--", "", "Output folder must be provided!!"),
    ],
)
def test_original_validation(sg, selections, monkeypatch, field, value, message):
    """Keep the original validation popups before invoking the launcher."""
    execute = Mock()
    monkeypatch.setattr("src.gui.run_nextflow", execute)
    selections[field] = value
    cancel_after(sg, [("--SUBMIT--", selections)])
    assert sg.PopupError.call_args.args == (message,)
    execute.assert_not_called()


def test_empty_fastq_folder(sg, selections, tmp_path):
    """Clear the input field after the original empty-folder error."""
    selections["--INDIR--"] = str(tmp_path)
    cancel_after(sg, [("--SUBMIT--", selections)])
    assert sg.PopupError.call_args.args == ("Input folder has no fastqs!!",)
    sg.Window.return_value.Element.assert_called_once_with("--INDIR--")


@pytest.mark.parametrize(
    "event, text",
    [
        ("--DOC--", "https://ccbr.github.io/RENEE_NF/"),
        ("--HELP--", "ccbr_pipeliner@mail.nih.gov"),
    ],
)
def test_clipboard_actions(sg, monkeypatch, event, text):
    """Keep the documentation/help clipboard and popup interactions."""
    clipboard = Mock()
    monkeypatch.setattr("src.gui.copy_to_clipboard", clipboard)
    cancel_after(sg, [(event, {})])
    clipboard.assert_called_once_with(text)
    assert text in sg.Popup.call_args.args[0]


def test_fastq_extensions(tmp_path):
    """Retain both original compressed FASTQ extensions."""
    for name in ("sample.R1.fastq.gz", "sample.R2.fq.gz", "ignore.txt"):
        (tmp_path / name).touch()
    assert {path.rsplit("/", 1)[-1] for path in get_fastqs(str(tmp_path))} == {
        "sample.R1.fastq.gz",
        "sample.R2.fq.gz",
    }


def test_preset_choices():
    """Populate the genome selector from Nextflow's bundled reference configs."""
    choices = genome_choices(repo_base())
    assert choices == sorted(set(choices))
    assert "hg38_36" in choices
    assert "hg38_stub" not in choices


def test_capture_launcher_output():
    """Capture combined stdout/stderr and the CLI exit status for the popups."""
    result = run_nextflow(
        [
            sys.executable,
            "-c",
            "import sys; print('preview'); print('failed', file=sys.stderr); sys.exit(7)",
        ]
    )
    assert result.returncode == 7
    assert "preview" in result.stdout
    assert "failed" in result.stdout


def test_command_works_outside_checkout(tmp_path, monkeypatch):
    """Resolve the checkout entry point without depending on the working directory."""
    monkeypatch.chdir(tmp_path)
    result = run_nextflow([sys.executable, str(repo_base("main.py")), "run", "--help"])
    assert result.returncode == 0
    assert "Run the workflow" in result.stdout


def test_launch_failure_is_visible():
    """Surface OS failures rather than presenting empty, successful output."""
    with pytest.raises(click.ClickException, match="Could not launch RENEE"):
        run_nextflow(["/nonexistent/renee-test-command"])


def test_gui_command_is_lazy(monkeypatch):
    """Only the GUI command loads its implementation."""
    launch = Mock()
    monkeypatch.setattr("src.gui.launch_gui", launch)
    result = CliRunner().invoke(cli, ["gui"])
    assert result.exit_code == 0, result.output
    launch.assert_called_once_with(repo_base())


def test_help_without_gui_dependencies():
    """Non-GUI commands must not import PySimpleGUI or Tkinter."""
    script = (
        "import sys; sys.modules['tkinter'] = None; "
        "sys.modules['PySimpleGUI'] = None; "
        "from src.__main__ import cli; cli(['gui', '--help'])"
    )
    subprocess.run(
        [sys.executable, "-c", script], cwd=repo_base(), check=True, capture_output=True
    )
