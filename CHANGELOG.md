## RENEE_NF development version

### API-breaking changes

- None.

### New features

- Implemented `RSeQC` for post-alignment RNA-seq quality control, including read distribution, infer experiment, and junction annotation analyses.
- Added a local `BAM2STRANDEDBW` module that converts STAR-aligned, duplicate-marked BAMs into forward and reverse strand BigWig files, porting the `bam2bw_rnaseq` rule from Snakemake RENEE. Handles both paired-end and single-end libraries and supports strand swapping via `--swap_strands` for non-dUTP libraries.
- STAR align now uses the max read length of samples to calculate --sjdbOverhang, matching snakemake
- Incorporated the nf-core `qualimap/bamqc` module to perform post-alignment BAM quality control. (#61, @AlecSilver)
- Added CCBR `samtools_flagstat` module into pipeline (#60, @AlecSilver)
- Added first steps of `picard` (add-or-replace-groups, and mark-duplicates) (#56, @AlecSilver)
- **RSEM integration**: Ported the `rsem` and `rsem_merge` rules from the RENEE Snakemake pipeline into Nextflow.
  - Installed `rsem/calculateexpression` from nf-core modules; patched container to `nciccbr/ccbr_rsem_1.3.3:v1.0` via `modules/nf-core/rsem/calculateexpression/rsem-calculateexpression.diff`.
  - Created local `RSEM_MERGE` module (`modules/local/rsem_merge/main.nf`) using the RENEE Python merge script (`bin/merge_rsem_results.py`), producing annotated expected-count, FPKM, TPM, and reformatted TSV matrices per the Snakemake version.
  - Created `subworkflows/local/rsem/main.nf` which: reads `strand.info` from `infer_experiment.py` to compute `meta.strandedness` (mirroring the Snakemake `--forward-prob` logic), stages the RSEM reference directory from the `rsem_ref` prefix in genome configs, calls `RSEM_CALCULATEEXPRESSION` in `--alignments` (BAM) mode, and gathers results into `RSEM_MERGE`.
  - Wired `RSEM` subworkflow into `main.nf` downstream of `STAR_ALIGN`; all outputs published to `DEG_ALL/`.
  - Parameters match Snakemake exactly: `--no-bam-output --calc-ci --seed 12345 --time --keep-intermediate-files --estimate-rspd`.
  - `rsem_ref` genome config key already present in all Biowulf and FRCE genome configs.
  - Added nf-test stub tests for `RSEM_CALCULATEEXPRESSION` and the `RSEM` subworkflow using sarscov2 data (both passing).
- **Agent & hooks**: Added `.github/agents/renee-nf-dev.agent.md` specialist agent for future RENEE_NF porting work; added `.github/hooks/nf-test-runner.json` + companion script that automatically runs `nf-test` with the Singularity profile after editing any `.nf.test` file.

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
