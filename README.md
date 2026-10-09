# RENEE_NF

**R**na s**E**quencing a**N**alysis pip**E**lin**E** in Nextflow

[![build](https://github.com/CCBR/RENEE_NF/actions/workflows/build-nextflow.yml/badge.svg)](https://github.com/CCBR/RENEE_NF/actions/workflows/build-nextflow.yml)
[![docs](https://github.com/CCBR/RENEE_NF/actions/workflows/docs-mkdocs.yml/badge.svg)](https://github.com/CCBR/RENEE_NF/actions/workflows/docs-mkdocs.yml)

> [!WARNING]
> This pipeline is under active development. We are accepting beta testers! Please report any issues you encounter (see [Help & Contributing](#help--contributing)).

This is a Nextflow port of [RENEE](https://ccbr.github.io/RENEE/latest/), an open-source, reproducible, and scalable pipeline for analyzing RNA-seq data.
See the website for detailed information, documentation, and examples:
<https://ccbr.github.io/RENEE_NF/>

## Introduction

RNA-sequencing (_RNA-seq_) has a wide variety of applications. This popular transcriptome profiling technique can be used to quantify gene and isoform expression, detect alternative splicing events, predict gene-fusions, call variants and much more.

**RENEE_NF** is a Nextflow port of [RENEE](https://github.com/CCBR/RENEE)[^1] (originally implemented in Snakemake[^2]), a comprehensive, open-source RNA-seq pipeline that relies on technologies like [Docker](https://www.docker.com/why-docker)[^3] and [Singularity/Apptainery](https://apptainer.org/docs/)[^4] to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Nextflow](https://docs.seqera.io/nextflow/)[^5], a flexible and scalable workflow management system, to submit jobs to a cluster or cloud provider.

## Overview

### RENEE Pipeline

A bioinformatics pipeline is more than the sum of its data processing steps. A pipeline without quality-control steps provides a myopic view of the potential sources of variation within your data (i.e., biological versus technical sources of variation). RENEE pipeline is composed of a series of quality-control and data processing steps.

The accuracy of the downstream interpretations made from transcriptomic data is highly dependent on the initial sample library. Unwanted sources of technical variation, if not accounted for properly, can influence the results. RENEE's comprehensive quality-control helps ensure your results are reliable and _reproducible across experiments_. In the data processing steps, RENEE quantifies gene and isoform expression and predicts gene fusions. Please note that the detection of alternative splicing events and variant calling will be incorporated in a later release.

![RNA-seq quantification pipeline](https://raw.githubusercontent.com/CCBR/RENEE/5e88c208846d8098a0e26570243cb987bc2d7402/resources/RENEE_Pipeline.svg) <sup>**Fig 1. An Overview of RENEE Pipeline.** Gene and isoform counts are quantified and a series of QC-checks are performed to assess the quality of the data. This pipeline stops at the generation of a raw counts matrix and gene-fusion calling. To run the pipeline, a user must select their raw data, a reference genome, and output directory (i.e., the location where the pipeline performs the analysis). Quality-control information is summarized across all samples in a MultiQC report.</sup>

**Quality Control**
[_FastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)[^6] is used to assess the sequencing quality. FastQC is run twice, before and after adapter trimming. It generates a set of basic statistics to identify problems that can arise during sequencing or library preparation. FastQC will summarize per base and per read QC metrics such as quality scores and GC content. It will also summarize the distribution of sequence lengths and will report the presence of adapter sequences.

[_Kraken2_](http://ccb.jhu.edu/software/kraken2/)[^7] and [_FastQ Screen_](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)[^8] are used to screen for various sources of contamination. During the process of sample collection to library preparation, there is a risk of introducing unwanted sources of DNA. FastQ Screen compares your sequencing data to a set of different reference genomes to determine if there is contamination. It allows a user to see if the composition of your library matches what you expect. Also, if there are high levels of microbial contamination, Kraken can provide an estimation of the taxonomic composition. Kraken can be used in conjunction with [_Krona_](https://github.com/marbl/Krona/wiki/KronaTools)[^9] to produce interactive reports.

[_Preseq_](http://smithlabresearch.org/software/preseq/)[^10] is used to estimate the complexity of a library for each sample. If the duplication rate is very high, the overall library complexity will be low. Low library complexity could signal an issue with library preparation where very little input RNA was over-amplified or the sample may be degraded.

[_Picard_](https://broadinstitute.github.io/picard/)[^11] can be used to estimate the duplication rate, and it has another particularly useful sub-command called CollectRNAseqMetrics which reports the number and percentage of reads that align to various regions: such as coding, intronic, UTR, intergenic and ribosomal regions. This is particularly useful as you would expect a library constructed with poly(A)-selection to have a high percentage of reads that map to coding regions. Picard CollectRNAseqMetrics will also report the uniformity of coverage across all genes, which is useful for determining whether a sample has a 3' bias (observed in poly(A)-selection libraries containing degraded RNA).

[_RSeQC_](http://rseqc.sourceforge.net/)[^12] is another particularly useful package that is tailored for RNA-seq data. It is used to calculate the inner distance between paired-end reads and calculate TIN values for a set of canonical protein-coding transcripts. A median TIN value is calculated for each sample, which is analogous to a computationally derived RIN.

[MultiQC](https://multiqc.info/)[^13] is used to aggregate the results of each tool into a single interactive report.

**Quantification**
[_Cutadapt_](https://cutadapt.readthedocs.io/en/stable/)[^14] is used to remove adapter sequences, perform quality trimming, and remove very short sequences that would otherwise multi-map all over the genome prior to alignment.

[_STAR_](https://github.com/alexdobin/STAR)[^15] is used to align reads to the reference genome. The RENEE pipeline runs STAR in two passes, where splice-junctions are collected and aggregated across all samples and provided to the second-pass of STAR. In the second pass of STAR, the splice-junctions detected in the first pass are inserted into the genome indices prior to alignment.

[_RSEM_](https://github.com/deweylab/RSEM)[^16] is used to quantify gene and isoform expression. The expected counts from RSEM are merged across samples to create two count matrices for gene counts and isoform counts.

[_Arriba_](https://arriba.readthedocs.io/en/latest/)[^17] is used to predict gene-fusion events. The pre-built human and mouse reference genomes use Arriba blacklists to reduce the false-positive rate.

Additional software, databases, and relevant citations:

- GENCODE[^18]
- samtools[^19]
- qualimap2[^20]
- R[^21]
- voom[^22]
- edgeR[^23]
- empirical Bayes[^24]
- fusion transcript detection benchmarking[^25]

### Reference Genomes

Pre-built reference genomes are provided on Biowulf and FRCE for a number of different annotation versions, view the list here:
<https://ccbr.github.io/RENEE/latest/RNA-seq/Resources/#1-reference-genomes>

If you would like to use a custom reference that is not already listed above,
you can prepare it with `renee_nf run --build`. See [Building a Custom Genome on Biowulf](#building-a-custom-genome-on-biowulf)
or [Setup for Generic SLURM Cluster](#setup-for-generic-slurm-cluster).

### Dependencies

**Requires:** `Nextflow>=25.10` and a container engine: `singularity>=3.5` / Apptainer or Docker

> [!NOTE]
>
> <ins>Biowulf users</ins>:
> Both the singularity and Nextflow modules are already installed and available for all Biowulf users. Please skip this step as `module load nextflow` will preload singularity and Nextflow.

[Nextflow](https://docs.seqera.io/nextflow/install) and a container engine, either [Singularity/Apptainer](https://apptainer.org/docs/admin/main/installation.html) or [Docker](https://docs.docker.com/get-docker/), must be installed on the target system. Nextflow orchestrates the execution of each step in the pipeline. To guarantee reproducibility, each step relies on pre-built images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). With the `singularity` profile, Nextflow pulls these Docker images, converts them to Singularity images on the fly, and saves them onto the local filesystem prior to job execution; with the `docker` profile, the images are run directly. As such, Nextflow and a container engine are the only two dependencies.

<hr>
<p align="center">
	<a href="#renee_nf">Back to Top</a>
</p>
<hr>

## Run the pipeline

### Samplesheet preparation

RENEE takes a comma-separated samplesheet, passed with `--input`, that lists every FastQ file to analyze. Each row is one sequencing library and the file must contain a header with the following columns:

| Column      | Required | Description                                                                                |
| ----------- | -------- | ------------------------------------------------------------------------------------------ |
| `replicate` | yes      | Replicate identifier (e.g. `S1`, `rep1`).                                                  |
| `sample`    | yes      | Sample or group name (e.g. `WT`, `KO`). Replicates of the same condition share this value. |
| `fastq_1`   | yes      | Path to the Read 1 FastQ file (typically gzipped, `.fastq.gz`).                            |
| `fastq_2`   | no       | Path to the Read 2 FastQ file for paired-end data. Leave empty for single-end data.        |

Example (paired-end):

```csv
replicate,sample,fastq_1,fastq_2
S1,WT,/path/to/WT_S1.R1.fastq.gz,/path/to/WT_S1.R2.fastq.gz
S2,WT,/path/to/WT_S2.R1.fastq.gz,/path/to/WT_S2.R2.fastq.gz
S3,KO,/path/to/KO_S3.R1.fastq.gz,/path/to/KO_S3.R2.fastq.gz
S4,KO,/path/to/KO_S4.R1.fastq.gz,/path/to/KO_S4.R2.fastq.gz
```

For single-end libraries, leave `fastq_2` empty but keep the trailing comma, e.g. `S4,KO,/path/to/KO_S4.R1.fastq.gz,`. Single-end and paired-end rows may be mixed in the same samplesheet.

Notes:

- Each library is named `<sample>_<replicate>` in the output (e.g. `WT_S1`), so every `sample`/`replicate` combination must be unique. Stick to letters, numbers, and underscores.
- Use absolute paths to FastQ files. `renee_nf run` launches Nextflow from the `--output` directory, so relative paths are resolved against that directory, not where you ran the command.
- All FastQ files are checked with fastQValidator before analysis; the run stops and lists any files that fail validation.
- Example samplesheets are available in [`assets/`](assets/) (`samplesheet.csv`, `samplesheet_single_read.csv`).

### Biowulf

Load the ccbrpipeliner module:

> [!NOTE]
>
> `renee_nf` will be available in `ccbrpipeliner` release 9 and later.

```sh
# grab an interactive node
sinteractive --mem=110g --cpus-per-task=12 --gres=lscratch:200
# load the module
module load ccbrpipeliner/9
```

RENEE supports two execution modes:

- slurm: submits each step as a job to the cluster (recommended)
- local: runs every step serially on the current node; useful for testing and debugging.

```sh
# View the help page for more information
renee_nf run --help
# Initialize the output directory
renee_nf init --output /data/$USER/RNA_hg38
cd /data/$USER/RNA_hg38
# Edit your samplesheet with your preferred text editor
nano assets/samplesheet.csv
# Launch the pipeline
renee_nf run \
    --mode slurm \
    --input assets/samplesheet.csv \
    --genome hg38_36 \
    --output /data/$USER/RNA_hg38

# Tip: add the -preview flag to see which processes will run
# without executing them.
renee_nf run \
    -preview \
    --mode local \
    --input assets/samplesheet.csv \
    --genome hg38_36 \
    --output /data/$USER/RNA_hg38
```

<!-- TODO: Add FRCE instructions -->

### Building a Custom Genome on Biowulf

If the genome or annotation you need is not one of the [pre-built reference genomes](#reference-genomes), you can build it once from a FASTA and GTF file and reuse it for every project. On Biowulf, the `biowulf` profile already points to the CCBR shared resources (Kraken2 and FastQ Screen databases) and SIF cache, so there is no need to pass `--download_shared_resources`.

```sh
sinteractive --mem=16g --cpus-per-task=2
module load ccbrpipeliner/9

renee_nf init \
    --output /data/$USER/renee_refs

renee_nf run \
    -profile biowulf,slurm \
    --build \
    --publish_dir_mode copy \
    --genome <genome name> \
    --genome_fasta /path/to/genome.fa \
    --genes_gtf /path/to/genes.gtf \
    --arriba_db_dir /data/CCBR_Pipeliner/db/PipeDB/arriba/arriba_v2.4.0/database \
    --output /data/$USER/renee_refs
```

When the build finishes, the genome resources and a `<genome name>.config` file are written to the output directory (`/data/$USER/renee_refs/genome/` in the example above).
Pass that config with `-c` to analyze your samples using the custom genome:

```sh
renee_nf init \
    --output /data/$USER/my_project

renee_nf run \
    -profile biowulf,slurm \
    --input samplesheet.csv \
    --genome <genome name> \
    -c /data/$USER/renee_refs/genome/<genome name>.config \
    --output /data/$USER/my_project
```

Gene-fusion calling with Arriba requires genome-specific reference files. If `<genome name>` contains `hg19`, `hg38`, `mm10`, or `mm39` (or an alias, e.g. `GRCh38`, `GRCm39`), they are detected automatically from `--arriba_db_dir`. For any other genome, either pass `--fusion_blacklist`, `--fusion_cytoband`, `--fusion_protdomain`, and `--fusion_known_fusions` explicitly, or omit them and fusion calling will be skipped.

### Setup for other platforms

Running the pipeline outside of Biowulf is easy; however, there are a few extra options you must provide.
Skip the below section if you are running the pipeline on Biowulf.

#### Installation on other platforms

You will need nextflow and a container engine such as docker, podman, or singularity/apptainer.

If you would like to use the `renee_nf` wrapper CLI, you will need python3 and then install it with pip or uv:

```sh
pip3 install git+https://github.com/CCBR/RENEE_NF
renee_nf --help
```

Alternatively, you can call the nextflow pipeline directly:

```sh
nextflow run CCBR/RENEE_NF -profile singularity,slurm ...
```

#### Build resources on other platforms

First, the genome must be built from the corresponding FASTA and GTF files. Additionally, when running the build (`--build`) for the first time, you will also need to provide the --download_shared_resources option. This option will download our kraken2 database and bowtie2 indices for FastQ Screen.

<!-- TODO: Mention setting up temp dir and SIF cache. -->

```sh
renee_nf init \
    --output path/to/ref_dir

renee_nf run \
  -profile <singularity | docker> \
  --mode <local | slurm> \
  --build \
  --publish_dir_mode copy \
  --download_shared_resources \
  --genome <genome name> \
  --genome_fasta /path/to/genome.fa \
  --genes_gtf /path/to/genome.gtf \
  --output path/to/ref_dir
```

If `<genome name>` (the value passed to `--genome`) contains `hg19`, `hg38`, `mm10`, or `mm39` (or an alias, e.g. `GRCh38`, `GRCm39`), RENEE automatically detects the matching Arriba fusion-calling reference files in the shared resources directory populated by `--download_shared_resources` above (or in the directory passed via `--arriba_db_dir`, if you set that instead). For any other genome, pass `--fusion_blacklist`, `--fusion_cytoband`, `--fusion_protdomain`, and `--fusion_known_fusions` explicitly.

### Building another genome

When building another genome, skip re-downloading the large shared resources files by omitting `--download_shared_resources` and instead including `-c path/to/ref_dir/shared_resources.config`.

```sh
renee_nf run \
    -profile <singularity | docker> \
    -c path/to/ref_dir/shared_resources.config \
    --mode <local | slurm> \
    --build \
    --publish_dir_mode copy \
    --genome <second genome name> \
    --genome_fasta /path/to/second_genome.fa \
    --genes_gtf /path/to/second_genome.gtf \
    --output path/to/ref_dir
```

### Running on Generic SLURM Cluster

Once the genome resources have been built, run the pipeline against your samples. Replace `path/to/ref_dir` with the output directory used to build the genome resources above, and `<genome name>` with the genome name chosen there.

```sh
renee_nf init \
  --output path/to/project

renee_nf run \
    -profile <singularity | docker> \
    --mode <local | slurm> \
    --input input.csv \
    --genome <genome name> \
    -c path/to/ref_dir/genome/<genome name>.config \
    -c path/to/ref_dir/shared_resources.config \
    --output path/to/project
```

The `-c` flag lets you include additional config files, each supplying extra parameters for the pipeline. Outside of Biowulf, this pipeline needs two additional configs:

- `<genome name>.config` adds paths specific to the genome being used.
- `shared_resources.config` points to the Kraken2, FastQ Screen, and Arriba databases, which are shared across all genome builds.

## About

This repo was originally generated from the [CCBR Nextflow Template](https://github.com/CCBR/CCBR_NextflowTemplate).
The template takes inspiration from nektool[^26] and the nf-core[^27] template, and relies on nf-core[^28] and custom modules[^29].

## Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/RENEE_NF/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/RENEE_NF/discussions) or contact [CCBR_Pipeliner@mail.nih.gov](mailto:CCBR_Pipeliner@mail.nih.gov).

Want to **contribute** to this project? Check out the [contributing guidelines](.github/contributing.md).

<hr>
<p align="center">
	<a href="#renee_nf">Back to Top</a>
</p>
<hr>

## References

[^1]: Sevilla, S., Sovacool, K., Kuhn, S., Tandon, M., Koparde, V. RENEE: Rna sEquencing aNalysis pipElinE (Snakemake version). <https://github.com/CCBR/RENEE>. DOI: <https://doi.org/10.5281/zenodo.10553198>.

[^2]: Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600. DOI: <https://doi.org/10.1093/bioinformatics/bty350>.

[^3]: Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2.

[^4]: Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459. DOI: <https://doi.org/10.1371/journal.pone.0177459>.

[^5]: Di Tommaso, P., et al. (2017). "Nextflow enables reproducible computational workflows." Nat Biotechnol 35(4): 316-319. DOI: <https://doi.org/10.1038/nbt.3820>.

[^6]: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.

[^7]: Wood, D. E. and S. L. Salzberg (2014). "Kraken: ultrafast metagenomic sequence classification using exact alignments." Genome Biol 15(3): R46. DOI: <https://doi.org/10.1186/gb-2014-15-3-r46>.

[^8]: Wingett, S. and S. Andrews (2018). "FastQ Screen: A tool for multi-genome mapping and quality control." F1000Research 7(2): 1338. DOI: <https://doi.org/10.12688/f1000research.15931.2>.

[^9]: Ondov, B. D., et al. (2011). "Interactive metagenomic visualization in a Web browser." BMC Bioinformatics 12(1): 385. DOI: <https://doi.org/10.1186/1471-2105-12-385>.

[^10]: Daley, T. and A.D. Smith, Predicting the molecular complexity of sequencing libraries. Nat Methods, 2013. 10(4): p. 325-7. DOI: <https://doi.org/10.1038/nmeth.2375>.

[^11]: The Picard toolkit. <https://broadinstitute.github.io/picard/>.

[^12]: Wang, L., et al. (2012). "RSeQC: quality control of RNA-seq experiments." Bioinformatics 28(16): 2184-2185. DOI: <https://doi.org/10.1093/bioinformatics/bts356>.

[^13]: Ewels, P., et al. (2016). "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32(19): 3047-3048. DOI: <https://doi.org/10.1093/bioinformatics/btw354>.

[^14]: Martin, M. (2011). "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet 17(1): 10-12. DOI: <https://doi.org/10.14806/ej.17.1.200>.

[^15]: Dobin, A., et al., STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 2013. 29(1): p. 15-21. DOI: <https://doi.org/10.1093/bioinformatics/bts635>.

[^16]: Li, B. and C.N. Dewey, RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics, 2011. 12: p. 323. DOI: <https://doi.org/10.1186/1471-2105-12-323>.

[^17]: Uhrig, S., et al. (2021). "Accurate and efficient detection of gene fusions from RNA sequencing data." Genome Res 31(3): 448-460. DOI: <https://doi.org/10.1101/gr.257246.119>.

[^18]: Harrow, J., et al., GENCODE: the reference human genome annotation for The ENCODE Project. Genome Res, 2012. 22(9): p. 1760-74. DOI: <https://doi.org/10.1101/gr.135350.111>.

[^19]: Li, H., et al. (2009). "The Sequence Alignment/Map format and SAMtools." Bioinformatics 25(16): 2078-2079. DOI: <https://doi.org/10.1093/bioinformatics/btp352>.

[^20]: Okonechnikov, K., et al. (2015). "Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data." Bioinformatics 32(2): 292-294. DOI: <https://doi.org/10.1093/bioinformatics/btv566>.

[^21]: R Core Team (2018). R: A Language and Environment for Statistical Computing. Vienna, Austria, R Foundation for Statistical Computing.

[^22]: Law, C.W., et al., voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol, 2014. 15(2): p. R29. DOI: <https://doi.org/10.1186/gb-2014-15-2-r29>.

[^23]: Robinson, M. D., et al. (2009). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics 26(1): 139-140. DOI: <https://doi.org/10.1093/bioinformatics/btp616>.

[^24]: Smyth, G.K., Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Stat Appl Genet Mol Biol, 2004. 3: p. Article3. DOI: <https://doi.org/10.2202/1544-6115.1027>.

[^25]: Haas, B. J., et al. (2019). "Accuracy assessment of fusion transcript detection via read-mapping and de novo fusion transcript assembly-based methods." Genome Biology 20(1): 213. DOI: <https://doi.org/10.1186/s13059-019-1842-9>.

[^26]: Roach, M. J., et al. (2022). "Ten simple rules and a template for creating workflows-as-applications." PLoS Computational Biology 18(12): e1010705. DOI: <https://doi.org/10.1371/journal.pcbi.1010705>. Nektool: <https://github.com/beardymcjohnface/nektool>.

[^27]: Ewels, P. A., et al. (2020). "The nf-core framework for community-curated bioinformatics pipelines." Nature Biotechnology 38(3): 276-278. DOI: <https://doi.org/10.1038/s41587-020-0439-x>.

[^28]: Langer, B. E., et al. (2025). "Empowering bioinformatics communities with Nextflow and nf-core." Genome Biology 26(1): 228. DOI: <https://doi.org/10.1186/s13059-025-03673-9>.

[^29]: Sovacool, K. and Koparde, V. (2026). "nf-modules: Reusable modules and subworkflows for CCBR Nextflow pipelines" (v0.2.0). Zenodo. DOI: <https://doi.org/10.5281/zenodo.10223357>. Repository: <https://github.com/CCBR/nf-modules>.
