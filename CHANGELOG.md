## RENEE_NF development version


### API-breaking changes

- None.

### New features

- Updated template placeholders and project metadata for RENEE.
- Added a samplesheet-derived paired-read input channel from `params.input` and pointed the bundled samplesheet at the test FASTQs.
- Installed the nf-core `cutadapt` module, tracked it in `modules.json`, and wired `CUTADAPT(ch_reads)` into the main workflow after the local `FASTQC` step.
- Added nested `snakemake_results` test-data links for the local RENEE hg38 Snakemake expected-output set.

### Bug fixes

- Normalized volatile Cutadapt CPU-count log lines in module snapshots so harmless `--cores` differences do not fail nf-test comparisons.

## RENEE_NF v0.1.0


This is the first release of RENEE_NF 🎉
