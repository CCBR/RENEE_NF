## RENEE_NF development version

### API-breaking changes

- None.

### New features

- Updated template placeholders and project metadata for RENEE.
- Added a samplesheet-derived paired-read input channel from `params.input` and pointed the bundled samplesheet at the test FASTQs.
- Installed the nf-core `cutadapt` module, tracked it in `modules.json`, and wired `CUTADAPT(ch_reads)` into the main workflow after the local `FASTQC` step.
- Added a local BBtools `BBTOOLS_BBMERGE` module for paired-end insert-size histograms and a matching nf-test module test covering Snakemake-compatible `*_insert_sizes.txt` output naming.
- Added nested `snakemake_results` test-data links for the local RENEE hg38 Snakemake expected-output set.
- Added a local `FASTQ_SCREEN` module and wired it into the main workflow (called twice on trimmed reads with separate configs for multi-organism and vector/rRNA screening), with HPC config files for Biowulf and FRCE and matching nf-tests.

### Bug fixes

- Normalized volatile Cutadapt CPU-count log lines in module snapshots so harmless `--cores` differences do not fail nf-test comparisons.
- Fixed local `FASTQC` module metadata scoping by switching to a closure-based `tag` directive and removing module-local `publishDir`, resolving `ERROR ~ No such variable: meta` during Nextflow preview/CI runs.

## RENEE_NF v0.1.0

This is the first release of RENEE_NF 🎉
