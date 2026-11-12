"""
Entrypoint for RENEE CLI

Check out the wiki for a detailed look at customizing this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import pathlib

import ccbr_tools.pipeline.util
import ccbr_tools.pkg_util
import cffconvert.cli.cli
import click

from .launcher import cache_images, execute, prepare_args


def repo_base(*paths):
    """Locate a file in the installed pipeline."""
    basedir = pathlib.Path(__file__).absolute().parent.parent
    return basedir.joinpath(*paths)


def print_citation_flag(ctx, param, value):
    """Print the citation for the eager global citation flag."""
    if not value or ctx.resilient_parsing:
        return
    ccbr_tools.pkg_util.print_citation(
        citation_file=repo_base("CITATION.cff"), output_format="bibtex"
    )
    ctx.exit()


@click.group(
    cls=ccbr_tools.pkg_util.CustomClickGroup,
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(
    ccbr_tools.pkg_util.get_version(repo_base=repo_base),
    "-v",
    "--version",
    is_flag=True,
)
@click.option(
    "-c",
    "--citation",
    callback=print_citation_flag,
    is_eager=True,
    is_flag=True,
    expose_value=False,
    help="Print the citation in bibtex format and exit.",
)
def cli():
    """RENEE: RNA-seq analysis pipeline in Nextflow.

    docs: https://ccbr.github.io/RENEE_NF

    For more options, run:
    renee [command] --help"""


help_msg_extra = """
\b
Nextflow options:
-profile <profile>    Nextflow profile to use (e.g. test)
-params-file <file>   Nextflow params file to use (e.g. assets/params.yml)
-preview              Preview the processes that will run without executing them
--input <csv|FASTQs>  Samplesheet CSV or sample.R1/R2.fastq.gz files
--genome <name>       Reference genome (default: hg38_36)
--small_rna           Small RNA analysis (single-end only)
--star_2_pass_basic   Per-sample STAR two-pass alignment
--shared_resources   Shared FastQ Screen/Kraken2/Arriba resources directory
--partition <name>   SLURM partition for driver and tasks
--time <walltime>    SLURM walltime for driver and tasks
--max_jobs <int>     Maximum queued SLURM tasks
--threads <int>      Total CPUs available in local mode
--tmp_dir <path>     Task scratch directory
--singularity_cache <path>  Private Singularity download cache
--sif_cache <path>   Nextflow SIF image cache
--create_nidap_folder  Copy portable reports and FPKM tables to NIDAP/

\b
EXAMPLES:
Execute with slurm:
    renee run --output path/to/outdir --mode slurm
Preview the processes that will run:
    renee run --output path/to/outdir --mode local -preview
Add nextflow args (anything supported by `nextflow run`):
    renee run --output path/to/outdir --mode slurm -profile test
    renee run --output path/to/outdir --mode slurm -profile test -params-file assets/params.yml
"""


@click.command(
    epilog=help_msg_extra,
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
    },
)
@click.option(
    "--main",
    "main_path",
    help="Path to the RENEE main.nf file or the GitHub repo (CCBR/RENEE_NF). Defaults to the version installed in the $PATH.",
    type=str,
    default=repo_base("main.nf"),
    show_default=True,
    hidden=True,
)
@click.option(
    "--output",
    help="Output directory path for renee init & run. Equivalient to nextflow launchDir. Defaults to your current working directory.",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=pathlib.Path.cwd(),
    show_default=False,
)
@click.option(
    "--mode",
    "_mode",
    help="Run mode (slurm, local)",
    type=click.Choice(["slurm", "local"]),
    default="slurm",
    show_default=True,
)
@click.option(
    "--forceall",
    "-F",
    "force_all",
    help="Force all processes to run (i.e. do not use nextflow -resume)",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.argument("nextflow_args", nargs=-1)
@click.option(
    "--dry-run", "--dry_run", is_flag=True, help="Preview without executing tasks."
)
@click.option("--runmode", type=click.Choice(["init", "run"]), default="run")
@click.option("--wait", is_flag=True, help="Wait for the SLURM driver job to complete.")
@click.option("--partition", help="SLURM partition for driver and tasks.")
@click.option("--time", "walltime", help="SLURM walltime for driver and tasks.")
@click.option(
    "--singularity-cache",
    "--singularity_cache",
    type=click.Path(file_okay=False),
    help="Private Singularity/Apptainer download cache.",
)
@click.pass_context
def run(
    ctx,
    main_path,
    output,
    _mode,
    force_all,
    nextflow_args,
    dry_run,
    runmode,
    wait,
    partition,
    walltime,
    singularity_cache,
):
    """
    Run the workflow

    New output directories are initialized automatically. Pipeline parameters use
    underscores instead of hyphens; see the Snakemake compatibility guide.

    docs: https://ccbr.github.io/RENEE_NF
    """
    if main_path != "CCBR/RENEE_NF" and not os.path.isfile(main_path):
        raise click.BadParameter(
            f"Path to the RENEE main.nf file not found: {main_path}"
        )
    output_dir = pathlib.Path(output).expanduser().resolve()
    if not (output_dir / "nextflow.config").is_file():
        ctx.invoke(init, output=output_dir)
    (output_dir / "log").mkdir(exist_ok=True)
    if runmode == "run":
        args = prepare_args(nextflow_args, output_dir, repo_base())
        if dry_run and "-preview" not in args:
            args.append("-preview")
        if partition:
            args.extend(["--partition", partition])
        if walltime:
            args.extend(["--time", walltime])
        execute(
            str(pathlib.Path(main_path).resolve())
            if main_path != "CCBR/RENEE_NF"
            else main_path,
            output_dir,
            _mode,
            force_all,
            args,
            wait,
            partition,
            walltime,
            singularity_cache,
        )


@click.command()
@click.option(
    "--output",
    help="Output directory path for renee init & run. Equivalient to nextflow launchDir. Defaults to your current working directory.",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=pathlib.Path.cwd(),
    show_default=False,
)
def init(output, **kwargs):
    """Initialize the working directory by copying the system default config files"""
    output_dir = output if isinstance(output, pathlib.Path) else pathlib.Path(output)
    ccbr_tools.pkg_util.msg_box(f"Initializing RENEE in {output_dir}")
    (output_dir / "log/").mkdir(parents=True, exist_ok=True)
    paths = ("nextflow.config", "conf/", "assets/")
    ccbr_tools.pipeline.util.copy_config(paths, repo_base=repo_base, outdir=output_dir)


@click.command(
    context_settings={
        "help_option_names": ["-h", "--help"],
        "ignore_unknown_options": True,
    },
    epilog=help_msg_extra,
)
@click.option("--output", type=click.Path(file_okay=False), default=".")
@click.option("--mode", type=click.Choice(["slurm", "local"]), default="slurm")
@click.option("--dry-run", "--dry_run", is_flag=True)
@click.option("--wait", is_flag=True)
@click.option("--partition")
@click.option("--time", "walltime")
@click.option(
    "--singularity-cache", "--singularity_cache", type=click.Path(file_okay=False)
)
@click.argument("nextflow_args", nargs=-1)
@click.pass_context
def build(
    ctx,
    output,
    mode,
    dry_run,
    wait,
    partition,
    walltime,
    singularity_cache,
    nextflow_args,
):
    """Build references using --ref_fa, --ref_gtf, --ref_name and --gtf_ver.

    These are native pipeline parameters. --shared_resources optionally downloads
    shared databases. Use run --build_shared_resources_only for databases alone.
    """
    ctx.invoke(
        run,
        output=output,
        _mode=mode,
        dry_run=dry_run,
        wait=wait,
        partition=partition,
        walltime=walltime,
        singularity_cache=singularity_cache,
        nextflow_args=("--build", *nextflow_args),
    )


@click.command()
@click.option(
    "--sif-cache", "--sif_cache", required=True, type=click.Path(file_okay=False)
)
@click.option("--dry-run", "--dry_run", is_flag=True)
def cache(sif_cache, dry_run):
    """Cache workflow containers in Nextflow's SIF cache (requires Nextflow and Singularity)."""
    cache_images(repo_base(), pathlib.Path(sif_cache).expanduser().resolve(), dry_run)


@click.command()
def debug():
    """Print the installed pipeline base directory and version."""
    click.echo(f"RENEE BASE: {repo_base()}")
    click.echo(ccbr_tools.pkg_util.get_version(repo_base=repo_base))


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
def gui():
    """Select inputs, preview, and submit a workflow using a desktop GUI."""
    from .gui import launch_gui

    launch_gui(repo_base())


@click.command()
@click.option("--output", required=True, type=click.Path(exists=True, file_okay=False))
def unlock(output):
    """Explain Nextflow's automatic lock handling; never remove a live session lock."""
    click.echo(
        f"Nextflow releases session locks when its driver exits; no unlock is needed for {output}. "
        "Stop any running driver before restarting. If a stale filesystem lock persists, "
        "check for live Java processes and contact your filesystem administrator."
    )


@click.command()
@click.argument(
    "citation_file",
    type=click.Path(exists=True),
    required=True,
    default=repo_base("CITATION.cff"),
)
@click.option(
    "--output-format",
    "-f",
    default="bibtex",
    help="Output format for the citation",
    type=cffconvert.cli.cli.options["outputformat"]["type"],
)
def cite(citation_file, output_format):
    """
    Print the citation in the desired format

    citation_file : Path to a file in Citation File Format (CFF) [default: the CFF for ccbr_tools]
    """
    ccbr_tools.pkg_util.print_citation(
        citation_file=citation_file, output_format=output_format
    )


cli.add_command(run)
cli.add_command(init)
cli.add_command(cite)
cli.add_command(build)
cli.add_command(cache)
cli.add_command(debug)
cli.add_command(unlock)
cli.add_command(gui)


def main():
    """Run the RENEE command-line interface."""
    cli(prog_name="renee")


if __name__ == "__main__":
    main()
