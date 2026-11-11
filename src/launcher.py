"""Launcher helpers for legacy FASTQ inputs and reliable Nextflow execution."""

import csv
import glob
import json
import os
import re
import subprocess
from pathlib import Path

import ccbr_tools.pipeline.nextflow
import click


def prepare_args(args: tuple[str, ...], output: Path, base: Path) -> list[str]:
    """Resolve CLI paths before changing launch directory; convert FASTQ lists to CSV."""
    schema = json.loads((base / "nextflow_schema.json").read_text())
    path_options = {"-c", "-config", "-params-file", "-work-dir", "-w"}
    for group in schema["$defs"].values():
        path_options.update(
            f"--{name}"
            for name, spec in group.get("properties", {}).items()
            if spec.get("format") in {"file-path", "directory-path"}
        )
    result = []
    index = 0
    while index < len(args):
        token = args[index]
        option, equals, value = token.partition("=")
        if option == "--input":
            inputs = [value] if equals else []
            index += 1
            while index < len(args) and not args[index].startswith("-"):
                inputs.append(args[index])
                index += 1
            if not inputs:
                raise click.UsageError("--input requires a samplesheet or FASTQ files.")
            if len(inputs) == 1 and inputs[0].endswith(".csv"):
                samplesheet = Path(inputs[0]).expanduser().resolve()
                if not samplesheet.is_file():
                    raise click.BadParameter(f"Samplesheet not found: {samplesheet}")
            else:
                files = []
                for pattern in inputs:
                    matches = sorted(glob.glob(os.path.expanduser(pattern)))
                    if not matches:
                        raise click.BadParameter(f"No FASTQ files match: {pattern}")
                    files.extend(Path(match).resolve() for match in matches)
                samplesheet = output / "log" / "launcher_samplesheet.csv"
                write_samplesheet(files, samplesheet)
            result.extend(["--input", str(samplesheet)])
        elif option == "--genome":
            if not equals:
                index += 1
                if index == len(args) or args[index].startswith("-"):
                    raise click.UsageError(
                        "--genome requires a name or reference config."
                    )
                value = args[index]
            if value.endswith((".json", ".config")):
                source = Path(value).expanduser().resolve()
                if not source.is_file():
                    raise click.BadParameter(f"Genome config not found: {source}")
                config = (
                    import_genome(source, output)
                    if source.suffix == ".json"
                    else source
                )
                result.extend(["-c", str(config), "--genome", source.stem])
            else:
                result.extend(["--genome", value])
            index += 1
        elif option in path_options:
            if not equals:
                index += 1
                if index == len(args) or args[index].startswith("-"):
                    raise click.UsageError(f"{option} requires a path.")
                value = args[index]
            result.extend([option, str(Path(value).expanduser().resolve())])
            index += 1
        else:
            result.append(token)
            index += 1
    return result


def import_genome(source: Path, output: Path) -> Path:
    """Translate a legacy references.rnaseq JSON file into a Nextflow genome preset."""
    keys = {
        "GENOMEFILE": "fasta",
        "GTFFILE": "genes_gtf",
        "GENOME_STARDIR": "star_index",
        "ANNOTATE": "annotate",
        "ANNOTATEISOFORMS": "annotate_isoforms",
        "REFFLAT": "refflat",
        "BEDREF": "bed_ref",
        "GENEINFO": "gene_info",
        "QUALIMAP_INFO": "qualimap_info",
        "KARYOBEDS": "karyobeds",
        "KARYOPLOTER": "karyoploter",
        "RSEMREF": "rsem_ref",
        "RRNALIST": "rrna_list",
        "ORGANISM": "organism",
        "TINREF": "tin_ref",
        "FUSIONBLACKLIST": "fusion_blacklist",
        "FUSIONCYTOBAND": "fusion_cytoband",
        "FUSIONPROTDOMAIN": "fusion_protdomain",
    }
    try:
        legacy = json.loads(source.read_text())["references"]["rnaseq"]
        missing = {
            "GENOMEFILE",
            "GTFFILE",
            "GENEINFO",
            "RSEMREF",
            "ORGANISM",
        } - legacy.keys()
        if missing:
            raise ValueError(
                f"Missing legacy reference fields: {', '.join(sorted(missing))}"
            )
        values = {new: legacy[old] for old, new in keys.items() if legacy.get(old)}
        if not all(isinstance(value, str) for value in values.values()):
            raise ValueError("Reference values must be strings.")
        values = {
            key: value
            if key == "organism"
            else str(Path(os.path.expandvars(value)).expanduser().resolve())
            for key, value in values.items()
        }
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise click.BadParameter(f"Cannot import genome {source}: {error}") from error

    def quote(value: str) -> str:
        """Quote a literal without Groovy interpolation."""
        return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"

    entries = ",\n".join(f"    {key}: {quote(value)}" for key, value in values.items())
    destination = output / "log" / "imported_genome.config"
    destination.write_text(f"params.genomes[{quote(source.stem)}] = [\n{entries}\n]\n")
    return destination


def write_samplesheet(files: list[Path], destination: Path) -> None:
    """Group legacy sample.R1/R2.fastq.gz files, rejecting duplicate or missing mates."""
    samples = {}
    seen = set()
    for fastq in files:
        match = re.fullmatch(r"(.+)\.R([12])\.(?:fastq|fq)\.gz", fastq.name)
        if not fastq.is_file() or not match:
            raise click.BadParameter(
                f"Expected an existing sample.R1.fastq.gz or sample.R2.fastq.gz: {fastq}"
            )
        if fastq in seen:
            raise click.BadParameter(f"Duplicate FASTQ input: {fastq}")
        seen.add(fastq)
        sample, mate = match.groups()
        reads = samples.setdefault(sample, {})
        if mate in reads:
            raise click.BadParameter(f"Multiple R{mate} files for sample {sample}")
        reads[mate] = str(fastq)
    if not samples:
        raise click.BadParameter("No FASTQ inputs provided.")
    paired = any("2" in reads for reads in samples.values())
    for sample, reads in samples.items():
        if "1" not in reads or (paired and "2" not in reads):
            raise click.BadParameter(f"Missing FASTQ mate for sample {sample}")
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "replicate", "fastq_1", "fastq_2"])
        for sample, reads in sorted(samples.items()):
            writer.writerow([sample, "", reads["1"], reads.get("2", "")])


def execute(
    main_path: str,
    output: Path,
    mode: str,
    force_all: bool,
    args: list[str],
    wait: bool = False,
    partition: str | None = None,
    walltime: str | None = None,
    singularity_cache: str | None = None,
) -> None:
    """Validate RENEE's private cache and delegate execution to the shared runner."""
    try:
        env = {}
        if singularity_cache:
            cache_dir = Path(singularity_cache).expanduser().resolve()
            cache_dir.mkdir(parents=True, exist_ok=True)
            if cache_dir.stat().st_uid != os.getuid():
                raise click.BadParameter(
                    "The private Singularity cache must be owned by you."
                )
            env = {
                "SINGULARITY_CACHEDIR": str(cache_dir),
                "APPTAINER_CACHEDIR": str(cache_dir),
            }
        ccbr_tools.pipeline.nextflow.run(
            nextfile_path=main_path,
            cwd=output,
            mode=mode,
            force_all=force_all,
            nextflow_args=args,
            pipeline_name="RENEE",
            hpc_memory="4G",
            hpc_walltime=walltime or "1-00:00:00",
            hpc_partition=partition,
            wait=wait,
            env=env,
            preview=False,
        )
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise click.ClickException(f"Workflow launch failed: {error}") from error


def cache_images(base: Path, destination: Path, dry_run: bool) -> None:
    """Inspect resolved Nextflow containers and populate its Singularity SIF cache."""
    try:
        inspection = subprocess.run(
            [
                "nextflow",
                "inspect",
                str(base / "main.nf"),
                "-profile",
                "singularity",
                "-format",
                "json",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        processes = json.loads(inspection.stdout)["processes"]
        images = sorted(
            {process["container"] for process in processes if process.get("container")}
        )
        if not images:
            raise click.ClickException(
                "Nextflow inspection did not report any containers."
            )
        if not dry_run:
            destination.mkdir(parents=True, exist_ok=True)
        for image in images:
            # Match Nextflow's SingularityCache naming convention.
            if image.startswith(("http://", "https://")):
                name = image.split("://", 1)[1].replace("/", "-")
                if not name.endswith(".sif"):
                    name += ".img"
                source = image
            else:
                source = image if image.startswith("docker://") else f"docker://{image}"
                name = (
                    source.removeprefix("docker://").replace("/", "-").replace(":", "-")
                    + ".img"
                )
            target = destination / name
            click.echo(f"{source} -> {target}")
            if not dry_run and not target.exists():
                subprocess.run(["singularity", "pull", str(target), source], check=True)
    except (OSError, subprocess.CalledProcessError, ValueError, KeyError) as error:
        detail = (
            error.stderr
            if isinstance(error, subprocess.CalledProcessError)
            else str(error)
        )
        raise click.ClickException(f"Container caching failed: {detail}") from error
