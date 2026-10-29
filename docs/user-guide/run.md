# 3. Running the Pipeline

RENEE 3 and later run with Nextflow. Use the installed `renee` command,
or `bin/renee_nf` from a checkout with the Python dependencies activated.
For every legacy option and the changes in its behavior, see
[Migrating from Snakemake](../snakemake-backward-compatibility.md).

## 3.1 Preview and run

```bash
# Preview without executing tasks or submitting a SLURM driver.
renee run --input reads/*.R?.fastq.gz --genome hg38_36 \
    --output analysis --dry-run

# Submit the driver and tasks to SLURM, waiting for completion.
renee run --input reads/*.R?.fastq.gz --genome hg38_36 \
    --output analysis --mode slurm --wait

# Or run locally with a container profile and explicit resource budgets.
renee run --input samplesheet.csv --genome hg38_36 \
    --output analysis --mode local --threads 4 --max_memory '12 GB' \
    -profile singularity
```

New launch directories are initialized automatically; `renee init --output analysis`
is optional. Results default to `analysis/results/`. Change the publication root
with `--outputDir`; change persistent task storage with `-work-dir`.
On a generic cluster or local system, supply a reference config if no genome
preset is available, for example `--genome reference-build/results/genome/hg38_47.config`.

## 3.2 Parameters and resuming

Pipeline parameters use double dashes and underscores: `--small_rna`,
`--star_2_pass_basic`, `--shared_resources`, and `--create_nidap_folder`.
Nextflow engine options use single dashes: `-profile`, `-params-file`, and `-c`.
See [pipeline parameters](params.md) for configuration-file usage.

Runs use Nextflow `-resume` automatically. Preserve both the `.nextflow/` session
metadata and the task work directory. Use `--forceall` (or `-F`) to disable
resume and execute tasks afresh. Do not reuse a Snakemake working directory.

## 3.3 Other commands

- `renee build`: build a reference from FASTA and GTF files.
- `renee cache`: populate a Nextflow-compatible Singularity image cache.
- `renee debug`: show the installed pipeline directory and version.
- `renee unlock --output analysis`: explain automatic Nextflow lock handling;
  no session files are deleted.
- `renee cite`: print the pipeline citation.

Use `renee COMMAND --help` for launcher options and
[the migration guide](../snakemake-backward-compatibility.md) for complete
reference-building, caching, and scheduler examples.
