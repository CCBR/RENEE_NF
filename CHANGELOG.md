## RENEE_NF development version

### API-breaking changes

- None.

### New features

- Added first steps of `picard` (add-or-replace-groups, and mark-duplicates) (#56, @AlecSilver)
- Added `PICARD_COLLECTRNASEQMETRICS` to the main workflow to collect RNA-seq alignment metrics (strand specificity, 5'/3' bias, UTR/intronic/coding distributions) from the duplicate-marked BAM, matching the RENEE Snakemake `stats` rule.
- Added `arriba` gene-fusion subworkflow (`subworkflows/local/arriba/main.nf`) with a dedicated `STAR_ALIGN_ARRIBA` modules.config block carrying all chimeric-detection flags (`--twopassMode Basic`, `--chimSegmentMin`, `--chimOutType WithinBAM HardClip`, etc.) matching the RENEE Snakemake `arriba` rule; updated nf-test module tests for `arriba/arriba` and `arriba/visualisation` to use local `params.test_data` references instead of `params.modules_testdata_base_path`.
- Added star_2_pass_basic mode, which does not pool splice junctions between samples in STAR alignment
- Updated template placeholders and project metadata for RENEE.
- Added a samplesheet-derived paired-read input channel from `params.input` and pointed the bundled samplesheet at the test FASTQs.
- Installed the nf-core `cutadapt` module, tracked it in `modules.json`, and wired `CUTADAPT(ch_reads)` into the main workflow after the local `FASTQC` step.
- Added a local BBtools `BBTOOLS_BBMERGE` module for paired-end insert-size histograms and a matching nf-test module test covering Snakemake-compatible `*_insert_sizes.txt` output naming.
- Incorporated the nf-core `fastq_screen` module to screen reads against multiple reference genomes and detect sample contamination.
- Added nested `snakemake_results` test-data links for the local RENEE hg38 Snakemake expected-output set.
- Migrated pipeline result publishing from process-level `publishDir` directives to entry-workflow `publish` and top-level `output` blocks, with dynamic per-sample paths and CSV index manifests for FastQC and Cutadapt outputs.

### Bug fixes

- Normalized volatile Cutadapt CPU-count log lines in module snapshots so harmless `--cores` differences do not fail nf-test comparisons.
- Fixed local `FASTQC` module metadata scoping by switching to a closure-based `tag` directive and removing module-local `publishDir`, resolving `ERROR ~ No such variable: meta` during Nextflow preview/CI runs.

## RENEE_NF v0.1.0

This is the first release of RENEE_NF 🎉
