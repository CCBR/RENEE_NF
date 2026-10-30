# Migrating from Snakemake RENEE

RENEE 3 and later use Nextflow; RENEE 2 and earlier use Snakemake.
The installed command is still `renee`. When running directly from a checkout,
use `bin/renee_nf` instead of the old checkout's `bin/renee`.
Activate an environment containing the Python dependencies before using either
the installed command or the checkout launcher. Nextflow >=25.10 and a container
engine are required for analysis; Snakemake is no longer required.

## Quick start

An old paired-end command:

```bash
renee run --input reads/*.R?.fastq.gz --output analysis \
    --genome hg38_36 --mode slurm --star-2-pass-basic --wait
```

becomes:

```bash
renee run --input reads/*.R?.fastq.gz --output analysis \
    --genome hg38_36 --mode slurm --star_2_pass_basic --wait
```

New working directories are initialized automatically. Existing configuration
files are not overwritten by `run`. You can still initialize explicitly with
`renee init --output analysis` or `renee run --output analysis --runmode init`.
Initialization only copies configuration; it does not save the other invocation
parameters. Supply analysis parameters when running, or save them in a Nextflow
params file and pass `-params-file params.yml`.

Pipeline parameters use **double dashes and underscores**.
Nextflow engine options use **single dashes**, such as `-profile`, `-resume`,
`-params-file`, `-c`, and `-work-dir`. Launcher options retain their existing
spelling, including `--output`, `--mode`, `--dry-run`, and `--wait`.
The launcher passes additional Nextflow options through without flattening
quoted arguments or discarding repeated configuration options.

## Complete option mapping

### Global options and commands

| Snakemake CLI         | Nextflow CLI                        | Notes                                                                                                                                                                                               |
| --------------------- | ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-h`, `--help`        | Unchanged                           | Also available for subcommands.                                                                                                                                                                     |
| `--version`           | Unchanged                           | `-v` is also available.                                                                                                                                                                             |
| `run`                 | Unchanged                           | Automatically initializes new working directories.                                                                                                                                                  |
| `build`               | Unchanged                           | Convenience command for `run --build`; reference flags are native parameters.                                                                                                                       |
| `cache`               | Unchanged                           | Uses Nextflow's container inspection and image naming; see cache differences below.                                                                                                                 |
| `debug`               | Unchanged                           | Prints installation directory and version.                                                                                                                                                          |
| `unlock --output DIR` | Same invocation, informational only | Nextflow releases locks when its driver exits. This command does not delete locks or session metadata. Stop live drivers before restarting; persistent filesystem locks require administrator help. |
| `gui`                 | Unchanged                           | Original PySimpleGUI form and dialogs, wired to Nextflow preview and SLURM submission; requires a graphical display.                                                                                |
| Not available         | `init`, `cite`, `--citation`        | Additional Nextflow launcher conveniences.                                                                                                                                                          |

### Analysis (`run`)

| Snakemake option          | Nextflow option                           | Behavior / differences                                                                                                                                                                                                                |
| ------------------------- | ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input FILE [FILE ...]` | Same, or `--input samplesheet.csv`        | Legacy FASTQ lists are converted to a samplesheet automatically. See input details below.                                                                                                                                             |
| `--output DIR`            | Unchanged                                 | Working/launch directory, not the published-results directory. Defaults to the current directory; see directory layout below.                                                                                                         |
| `--genome NAME`           | Unchanged                                 | Default is `hg38_36`; available presets depend on the profile. A legacy JSON or a generated Nextflow genome config path is also accepted by the launcher.                                                                             |
| `--small-rna`             | `--small_rna`                             | Single-end small RNA analysis. Cannot be combined with `--star_2_pass_basic`.                                                                                                                                                         |
| `--star-2-pass-basic`     | `--star_2_pass_basic`                     | Per-sample STAR two-pass mode rather than pooled splice junctions.                                                                                                                                                                    |
| `--dry-run`               | Same; `--dry_run` also accepted           | Maps to Nextflow `-preview`, without submitting a SLURM driver. It previews the workflow structure, not Snakemake's exact list of missing files/jobs.                                                                                 |
| `--runmode {init,run}`    | Unchanged                                 | `init` only initializes; `run` executes. Default: `run`.                                                                                                                                                                              |
| `--mode {slurm,local}`    | Unchanged                                 | Default: `slurm`. SLURM submits a driver with `sbatch`; local runs the driver in the foreground. HPC profiles are detected on supported clusters.                                                                                     |
| `--shared-resources DIR`  | `--shared_resources DIR`                  | Shared database root, reusable across output directories; explicit per-tool database/config parameters take precedence.                                                                                                               |
| `--singularity-cache DIR` | Same; `--singularity_cache` also accepted | Private download cache. Sets `SINGULARITY_CACHEDIR` and `APPTAINER_CACHEDIR`; the directory must be owned by you.                                                                                                                     |
| `--sif-cache DIR`         | `--sif_cache DIR`                         | Nextflow image cache; supports Singularity and Apptainer. Image names differ from the old Snakemake cache.                                                                                                                            |
| `--create-nidap-folder`   | `--create_nidap_folder`                   | Copies the RNA report, MultiQC report, gene FPKM table, and isoform FPKM table into `NIDAP/` under the published-results directory. Copies are portable and do not depend on work-directory links.                                    |
| `--wait`                  | Unchanged                                 | Adds `sbatch --wait` so the launcher waits for the SLURM driver and reports its failure. Local execution always waits.                                                                                                                |
| `--tmp-dir DIR`           | `--tmp_dir DIR`                           | Task scratch directory. Persistent Nextflow task/cache storage is controlled separately with `-work-dir DIR`; do not delete it if you want to resume.                                                                                 |
| `--partition NAME`        | Unchanged                                 | Applies to both the SLURM driver and tasks. Profile defaults apply when omitted.                                                                                                                                                      |
| `--time TIME`             | Unchanged                                 | SLURM walltime for the driver and maximum task walltime. Use `D-HH:MM:SS` or `HH:MM:SS` with the launcher. Direct Nextflow use also accepts durations such as `12 h` for tasks. Existing resource limits can further constrain tasks. |
| `--max-jobs N`            | `--max_jobs N`                            | Positive integer executor queue size; no automatic Snakemake student-partition 10-job default.                                                                                                                                        |
| `--threads N`             | Unchanged                                 | Positive integer local CPU budget and local per-task CPU cap; does not limit SLURM tasks. This replaces Snakemake's local `--cores` budget. Nextflow profile defaults apply when omitted.                                             |
| Not available             | `--forceall`, `-F`                        | Start without Nextflow `-resume`, forcing new task execution.                                                                                                                                                                         |

### Desktop launcher (`gui`)

Run `renee gui` or `bin/renee_nf gui` from a checkout. The original PySimpleGUI
implementation is retained, using the same `PySimpleGUI-4-foss < 5` dependency
as Snakemake RENEE. Python must have Tk support and a graphical display must be
available (for example, an X-forwarded session on Biowulf). CLI-only commands do
not import PySimpleGUI.

1. Select an input FASTQ folder using legacy `.R1`/`.R2` naming.
2. Select an output folder and a genome/annotation preset. The dropdown now reads
   the bundled Nextflow Biowulf/FRCE configs instead of Snakemake reference JSON
   files. References require their matching site profile and accessible files;
   supported HPC profiles are detected by the launcher.
3. Click **Submit** and review the scrollable preview output.
4. Confirm **Submit run to slurm?** to submit, or choose **No** to clear the
   selections and return to the form.

The original form, synchronous event loop, popup dialogs, clipboard actions, and
cancel exit status (69) are retained. The execution adapter invokes `renee run`
with `--mode slurm --threads 2`; the preview adds `--dry-run`, which uses
Nextflow `-preview`, not Snakemake's job-level dry run. Launcher exit statuses
replace the old error-text check because successful Nextflow previews include
process names such as `HANDLE_PRESEQ_ERROR`. Failed previews prevent submission;
submission failures are reported in an error popup.
After successful submission, inspect `log/slurm-*.out` under the launch directory
for workflow progress. Closing the GUI does not cancel submitted SLURM jobs.
**Documentation** now copies the Nextflow documentation URL; **Help** still
copies the support email address.

As in the original GUI, samplesheet input, custom reference browsing, local
execution, and additional analysis/resource controls are not exposed in the
form. Use `renee run` directly for those options.

### Reference building (`build`)

All orchestration flags in the legacy `build` command remain available:
`--output`, `--dry-run`, `--wait`, `--partition`, `--time`,
`--singularity-cache`, plus underscore versions of `--sif-cache` and `--tmp-dir`.
`--mode local` is additionally supported.

| Snakemake option         | Native Nextflow parameter | Notes                                                                                                                                                                        |
| ------------------------ | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--ref-fa FASTA`         | `--ref_fa FASTA`          | Alias for `--genome_fasta`.                                                                                                                                                  |
| `--ref-gtf GTF`          | `--ref_gtf GTF`           | Alias for `--genes_gtf`.                                                                                                                                                     |
| `--ref-name NAME`        | `--ref_name NAME`         | Combined with annotation version to produce `NAME_VERSION`.                                                                                                                  |
| `--gtf-ver VERSION`      | `--gtf_ver VERSION`       | Used in the built genome name.                                                                                                                                               |
| `--small-genome`         | `--small_genome`          | Accepted for compatibility, but no longer necessary: the nf-core STAR module automatically chooses a genome-size-dependent `genomeSAindexNbases` for every custom reference. |
| `--shared-resources DIR` | `--shared_resources DIR`  | In build mode, also enables shared database downloads to this directory. In analysis mode, reuses those databases.                                                           |

```bash
renee build --ref_fa references/genome.fa \
    --ref_gtf references/genes.gtf --ref_name hg38 --gtf_ver 47 \
    --output reference-build --mode slurm --wait

renee run --input reads/*.R?.fastq.gz --output analysis \
    --genome reference-build/results/genome/hg38_47.config --mode slurm
```

The equivalent native build invocation is `renee run --build` with the same
parameters. Advanced `--genome_fasta`/`--genes_gtf` spelling remains supported.
If both a native parameter and its legacy alias are given, the native value
takes precedence. Built reference descriptions are now Nextflow `.config`
files, not Snakemake `.json` files.

To download databases without building a reference:

```bash
renee run --build_shared_resources_only --shared_resources shared-databases \
    --output database-build --mode slurm --wait
```

The shared root contains `fastq_screen_db/`, `20180907_standard_kraken2/`,
and `arriba_v2.5.0/database/`. Builds also write `shared_resources.config`
under the published-results directory. Old shared resources can be reused
when their layout matches; otherwise supply `--fastq_screen_conf`,
`--fastq_screen_conf2`, `--fastq_screen_db_dir`, `--kraken2_db_dir`, and/or
`--arriba_db_dir` explicitly.

### Container caching (`cache`)

```bash
renee cache --sif-cache sif-cache --dry-run
renee cache --sif-cache sif-cache
renee run --input samplesheet.csv --output analysis --sif_cache sif-cache \
    --mode local -profile singularity
```

`--sif_cache` and `--dry_run` are also accepted by `cache`.
Dry-run inspects containers but does not pull images or create the cache directory.
Caching requires Nextflow; pulling additionally requires Singularity.
The cache uses Nextflow's image names, not the old Snakemake cache filenames,
so existing old cached images are not automatically reused.
This shared SIF cache is distinct from the user-owned private download cache.
Container inspection uses the Singularity profile.

## Input and sample identity

Legacy inputs must use `sample.R1.fastq.gz` / `sample.R2.fastq.gz` naming
(`.fq.gz` is also accepted). Expanded shell globs and quoted glob patterns both
work. Single-end runs supply only R1 files. Duplicates, missing mates, and
mixed single-/paired-end FASTQ lists are rejected before launch, matching the
legacy list-input restrictions.

The launcher writes `log/launcher_samplesheet.csv` with absolute FASTQ paths.
Each legacy basename becomes `sample`, with an empty `replicate`, preserving
the basename in result names. Explicit samplesheet rows with a replicate use
the existing Nextflow identity `sample_replicate`. Use a samplesheet for explicit sample and replicate metadata,
arbitrary FASTQ filenames, or mixed single-/paired-end experiments:

```csv
sample,replicate,fastq_1,fastq_2
tumor,1,/path/to/tumor_R1.fastq.gz,/path/to/tumor_R2.fastq.gz
control,1,/path/to/control.fastq.gz,
```

CLI file/directory paths are resolved relative to the invocation directory
before the launcher enters `--output`. Paths _inside_ a user-supplied CSV or
params file retain Nextflow's semantics: use absolute paths for portability;
relative FASTQ paths in a CSV are relative to the launch directory.
Bundled samplesheets use `${projectDir}` placeholders.

`--genome custom.json` imports the legacy `references.rnaseq` structure into
`log/imported_genome.config`. The original file is not changed, reference paths
must still be accessible, and container/tool-version compatibility is not
guaranteed merely because the JSON can be imported. A generated Nextflow
`--genome NAME.config` is shorthand for `-c NAME.config --genome NAME`.

## Output, profiles, and resuming

- `--output analysis` selects the launch directory. Published results default to
  `analysis/results/`; `--outputDir` controls this separately. For example,
  `--output analysis --outputDir .` publishes directly under `analysis`.
- The Nextflow `work/` directory, `.nextflow/` session database, trace files, and
  output layout replace Snakemake's working state. Start a **new output
  directory** when migrating; old Snakemake intermediates cannot be resumed.
- Runs resume automatically using Nextflow `-resume`. Keep both the session
  metadata and task work directories. `--forceall` disables resume.
- `--mode local` controls where the driver runs. To run tasks locally as well,
  do not select an executor profile such as `-profile slurm`. Container profiles
  such as `-profile docker` or `-profile singularity` are independent.
- Supported HPC profiles are added automatically. Generic SLURM clusters can
  use `--mode slurm`; configure storage, containers, and scheduler settings for
  the site. The launcher no longer loads Snakemake.
- Local workflow failures, submission failures, and waited-for SLURM driver
  failures return a nonzero launcher status. Without `--wait`, successful
  submission is not confirmation that the workflow completed successfully.

For all current pipeline parameters, see the
[parameter reference](user-guide/params.md) and `nextflow_schema.json`.
