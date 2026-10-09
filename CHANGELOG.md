## RENEE_NF development version

### API-breaking changes

- None.

### New features

- Build mode now downloads and creates all files for config, such as arriba files, and qualimap info and rsem refs
- Added `kraken2` module for taxonomic classification of reads. (#66, @AlecSilver)
- Added `MultiQC` report (#68, @AlecSilver)
- Implemented `RSeQC` for post-alignment RNA-seq quality control, including read distribution, infer experiment, and junction annotation analyses.
- Added a local `BAM2STRANDEDBW` module that converts STAR-aligned, duplicate-marked BAMs into forward and reverse strand BigWig files, porting the `bam2bw_rnaseq` rule from Snakemake RENEE. Handles both paired-end and single-end libraries and supports strand swapping via `--swap_strands` for non-dUTP libraries.
- Added `preseq` module from nf-core and incorporated into QC (#62, @AlecSilver)
- STAR align now uses the max read length of samples to calculate --sjdbOverhang, matching snakemake
- Incorporated the nf-core `qualimap/bamqc` module to perform post-alignment BAM quality control. (#61, @AlecSilver)
- Added CCBR `samtools_flagstat` module into pipeline (#60, @AlecSilver)
- Added first steps of `picard` (add-or-replace-groups, and mark-duplicates) (#56, @AlecSilver)
- Ported the `rsem` and `rsem_merge` rules from the RENEE Snakemake pipeline into Nextflow.
- Added `PICARD_COLLECTRNASEQMETRICS` to the main workflow to collect RNA-seq alignment metrics (strand specificity, 5'/3' bias, UTR/intronic/coding distributions) from the duplicate-marked BAM, matching the RENEE Snakemake `stats` rule. (#65, @AlecSilver)
- Added `arriba` gene-fusion subworkflow (`subworkflows/local/arriba/main.nf`) with a dedicated `STAR_ALIGN_ARRIBA` modules.config block carrying all chimeric-detection flags (`--twopassMode Basic`, `--chimSegmentMin`, `--chimOutType WithinBAM HardClip`, etc.) matching the RENEE Snakemake `arriba` rule; updated nf-test module tests for `arriba/arriba` and `arriba/visualisation` to use local `params.test_data` references instead of `params.modules_testdata_base_path`.
- Added `small_rna` mode, which trims reads to a minimum length of 16 bp and aligns them with STAR in a single pass followed by samtools sort, matching the `trim_se` and `star_small` rules in RENEE. (#92, @kelly-sovacool)
- Added star_2_pass_basic mode, which does not pool splice junctions between samples in STAR alignment
- Updated template placeholders and project metadata for RENEE.
- Added a samplesheet-derived paired-read input channel from `params.input` and pointed the bundled samplesheet at the test FASTQs.
- Installed the nf-core `cutadapt` module, tracked it in `modules.json`, and wired `CUTADAPT(ch_reads)` into the main workflow after the local `FASTQC` step.
- Added a local BBtools `BBTOOLS_BBMERGE` module for paired-end insert-size histograms and a matching nf-test module test covering Snakemake-compatible `*_insert_sizes.txt` output naming.
- Incorporated the nf-core `fastq_screen` module to screen reads against multiple reference genomes and detect sample contamination.
- Added nested `snakemake_results` test-data links for the local RENEE hg38 Snakemake expected-output set.
- Migrated pipeline result publishing from process-level `publishDir` directives to entry-workflow `publish` and top-level `output` blocks, with dynamic per-sample paths and CSV index manifests for FastQC and Cutadapt outputs.

### Bug fixes

- Fixed bundled samplesheet FASTQ paths to resolve relative to the pipeline directory. (#95, @kelly-sovacool)
- Fix `FASTQ_SCREEN_2` results being dropped by MultiQC, which crashed `RNA_REPORT`. (#91, @kelly-sovacool)
- Deduplicate `rNA_flowcells.Rmd` to match Snakemake RENEE. (#91, @kelly-sovacool)
- Fixed `bin/get_isoform_annotate.py` to derive isoform annotations from any GTF line carrying a `transcript_id` (CDS/exon/etc.), not just explicit `transcript`-type lines -- `gtfToGenePred` builds gene models from CDS/exon lines regardless, so a transcript_id with no matching `transcript` line (seen both with NCBI/GenBank-style GTFs and with hand-trimmed test fixtures that drop a transcript's header line) was silently missing from `annotate.isoforms.txt`, causing `make_refFlat.py` to crash with a `KeyError` in `BUILD_ANNOTATE`.
- Swapped container in Qualimap bamQC so behavior matches snakemake
- Normalized volatile Cutadapt CPU-count log lines in module snapshots so harmless `--cores` differences do not fail nf-test comparisons.
- Fixed local `FASTQC` module metadata scoping by switching to a closure-based `tag` directive and removing module-local `publishDir`, resolving `ERROR ~ No such variable: meta` during Nextflow preview/CI runs.

## RENEE_NF v0.1.0

This is the first release of RENEE_NF 🎉
