#!/usr/bin/env python3
"""Original RENEE PySimpleGUI launcher, wired to Nextflow.

Example: ``renee gui`` (requires a graphical display).
"""

import glob
import os
from pathlib import Path
import re
import subprocess
import sys

import click


def genome_choices(base: Path) -> list[str]:
    """Read the genome names used by the Nextflow reference configs."""
    names = set()
    for site in ("biowulf", "frce"):
        config = base / "conf" / "genomes" / f"{site}_genomes.config"
        names.update(
            re.findall(r"^\s*'([^']+)'\s*\{", config.read_text(), re.MULTILINE)
        )
    return sorted(names)


def launch_gui(base: Path, debug: bool = True) -> None:
    """Launch the legacy form and preview/SLURM submission dialogs."""
    # Deferred import: PySimpleGUI opens a network socket at import time.
    import PySimpleGUI as sg

    genome_annotation_combinations = genome_choices(base)
    if debug:
        print(genome_annotation_combinations)

    logo = sg.Image(str(base / "assets" / "CCBRlogo.png"))
    layout = [
        [sg.Column([[logo]], justification="center")],
        [
            sg.Text(
                "RENEE - Rna sEquencing aNalysis pipElinE", font=("Arial", 14, "bold")
            )
        ],
        [
            sg.Text(
                "Input Fastqs folder", font=("Helvetica", 12, "bold"), size=(20, 1)
            ),
            sg.InputText(key="--INDIR--"),
            sg.FolderBrowse(target="--INDIR--"),
        ],
        [
            sg.Text("Output folder", font=("Helvetica", 12, "bold"), size=(20, 1)),
            sg.InputText(key="--OUTDIR--"),
            sg.FolderBrowse(target="--OUTDIR--"),
        ],
        [
            sg.Text("Genome+Annotation", font=("Helvetica", 12, "bold"), size=(20, 1)),
            sg.Combo(
                values=genome_annotation_combinations,
                key="--ANNOTATION--",
                tooltip="eg. hg38_30 for Genome=hg38 & Gencode_Annotation=version 30",
            ),
        ],
        [
            sg.Submit(key="--SUBMIT--", font=("Helvetica", 12)),
            sg.Cancel(key="--CANCEL--", font=("Helvetica", 12)),
            sg.Button(
                button_text="Documentation", key="--DOC--", font=("Helvetica", 12)
            ),
            sg.Button(button_text="Help", key="--HELP--", font=("Helvetica", 12)),
        ],
    ]
    if debug:
        print("layout is ready!")

    window = sg.Window(
        f"RENEE {(base / 'VERSION').read_text().strip()}",
        layout,
        location=(0, 500),
        finalize=True,
    )
    if debug:
        print("window created!")

    while True:
        event, values = window.read()
        if debug:
            print(event, values)
        if event == "--CANCEL--" or event == sg.WIN_CLOSED:
            sg.popup_auto_close(
                "Thank you for running RENEE. GoodBye!",
                location=(0, 500),
                title="",
                font=("Arial", 12, "bold"),
            )
            sys.exit(69)
        if event == "--DOC--":
            copy_to_clipboard("https://ccbr.github.io/RENEE_NF/")
            sg.Popup(
                "Visit https://ccbr.github.io/RENEE_NF/ for links to complete documentation. The link has been copied to your clipboard. Please paste it in your favorite web browser.",
                font=("Arial", 12, "bold"),
                location=(0, 500),
            )
            continue
        if event == "--HELP--":
            copy_to_clipboard("ccbr_pipeliner@mail.nih.gov")
            sg.Popup(
                "Email ccbr_pipeliner@mail.nih.gov for help. The email id has been copied to your clipboard. Please paste it in your emailing software.",
                font=("Arial", 12, "bold"),
                location=(0, 500),
            )
            continue
        if event == "--SUBMIT--":
            if values["--INDIR--"] == "":
                sg.PopupError(
                    "Input folder must be provided!!",
                    location=(0, 500),
                    title="ERROR!",
                    font=("Arial", 12, "bold"),
                )
                continue
            elif not os.path.exists(values["--INDIR--"]) and not os.path.exists(
                fixpath(values["--INDIR--"])
            ):
                if debug:
                    print(values["--INDIR--"])
                if debug:
                    print(fixpath(values["--INDIR--"]))
                sg.PopupError(
                    "Input folder doesn't exist!!",
                    location=(0, 500),
                    title="ERROR!",
                    font=("Arial", 12, "bold"),
                )
                continue
            else:
                inputfastqs = get_fastqs(values["--INDIR--"])
                if debug:
                    print(inputfastqs)
                if len(inputfastqs) == 0:
                    sg.PopupError(
                        "Input folder has no fastqs!!",
                        location=(0, 500),
                        title="ERROR!",
                        font=("Arial", 12, "bold"),
                    )
                    window.Element("--INDIR--").update("")
                    continue
            if values["--OUTDIR--"] == "":
                sg.PopupError(
                    "Output folder must be provided!!",
                    location=(0, 500),
                    title="ERROR",
                    font=("Arial", 12, "bold"),
                )
                continue
            elif os.path.exists(values["--OUTDIR--"]) and not os.path.exists(
                fixpath(values["--OUTDIR--"])
            ):
                ch = sg.popup_yes_no(
                    "Output folder exists... this is probably a re-run ... proceed?",
                    title="Rerun?",
                    location=(0, 500),
                    font=("Arial", 12, "bold"),
                )
                if ch == "No":
                    window.Element("--OUTDIR--").update("")
                    continue
            genome = values["--ANNOTATION--"]
            output_dir = values["--OUTDIR--"]
            run_args = [
                sys.executable,
                str(base / "main.py"),
                "run",
                "--input",
                *inputfastqs,
                "--output",
                fixpath(output_dir),
                "--genome",
                genome,
                "--mode",
                "slurm",
                "--threads",
                "2",
                "--dry-run",
            ]
            result = run_nextflow(run_args)
            allout = result.stdout
            sg.popup_scrolled(
                allout,
                title="Dryrun:STDOUT/STDERR",
                font=("Monaco", 10),
                location=(0, 500),
                size=(80, 30),
            )
            # Nextflow previews include process names such as HANDLE_PRESEQ_ERROR.
            if result.returncode != 0:
                continue
            ch = sg.popup_yes_no(
                "Submit run to slurm?",
                title="Submit??",
                location=(0, 500),
                font=("Arial", 12, "bold"),
            )
            if ch == "Yes":
                result = run_nextflow(run_args[:-1])
                allout = result.stdout
                sg.popup_scrolled(
                    allout,
                    title="Dryrun:STDOUT/STDERR",
                    font=("Monaco", 10),
                    location=(0, 500),
                    size=(80, 30),
                )
                sg.popup_scrolled(
                    allout,
                    title="Slurmrun:STDOUT/STDERR",
                    font=("Monaco", 10),
                    location=(0, 500),
                    size=(80, 30),
                )
                if result.returncode != 0:
                    sg.PopupError(
                        f"Workflow submission failed (exit status {result.returncode}).",
                        location=(0, 500),
                        title="ERROR!",
                        font=("Arial", 12, "bold"),
                    )
                    continue
                sg.popup_auto_close(
                    "Thank you for running RENEE. GoodBye!",
                    location=(0, 500),
                    title="",
                    font=("Arial", 12, "bold"),
                )
                break
            elif ch == "No":
                window.Element("--INDIR--").update("")
                window.Element("--OUTDIR--").update("")
                window.Element("--ANNOTATION--").update("")
                continue

    window.close()


def run_nextflow(command: list[str]) -> subprocess.CompletedProcess[str]:
    """Capture the Nextflow launcher's output for the original GUI dialogs."""
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    except OSError as error:
        raise click.ClickException(f"Could not launch RENEE: {error}") from error
    return result


def copy_to_clipboard(string: str) -> None:
    """Copy the documentation URL or support email to the clipboard."""
    from tkinter import Tk

    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(string)
    r.update()
    r.destroy()


def fixpath(p: str) -> str:
    """Expand a GUI path to an absolute filesystem path."""
    return os.path.abspath(os.path.expanduser(p))


def get_fastqs(inputdir: str) -> list[str]:
    """Collect the same compressed FASTQ filenames as the Snakemake GUI."""
    inputdir = fixpath(inputdir)
    inputfastqs = glob.glob(inputdir + os.sep + "*.fastq.gz")
    inputfqs = glob.glob(inputdir + os.sep + "*.fq.gz")
    inputfastqs.extend(inputfqs)
    return inputfastqs
