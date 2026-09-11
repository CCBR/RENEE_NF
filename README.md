# RENEE_NF

bulk RNA-seq pipeline in Nextflow

[![build](https://github.com/CCBR/CCBR_NextflowTemplate/actions/workflows/build-nextflow.yml/badge.svg)](https://github.com/CCBR/CCBR_NextflowTemplate/actions/workflows/build-nextflow.yml)
[![docs](https://github.com/CCBR/CCBR_NextflowTemplate/actions/workflows/docs-mkdocs.yml/badge.svg)](https://github.com/CCBR/CCBR_NextflowTemplate/actions/workflows/docs-mkdocs.yml)

See the website for detailed information, documentation, and examples:
<https://ccbr.github.io/RENEE_NF/>

## Usage

Install the tool in edit mode:

```sh
pip3 install -e .
```

View CLI options:

```sh
renee_nf --help
```

Navigate to your project directory and initialize required config files:

```sh
renee_nf init
```

Run the example

```sh
renee_nf run --input "Hello world"
```

## Running on a new system

If you are running on a system other than Biowulf or FRCE, you will first need to build the genome resources and download the reference databases the pipeline requires.

### Build the genome resources

Initialize the output directory, then run the build:

```
renee_nf run \
    -profile <singularity | docker> \
    --mode <local | slurm> \
    --build \
    --publish_dir_mode copy \
    --genome <genome name> \
    --genome_fasta /path/to/genome.fa \
    --genes_gtf /path/to/.gtf \
    --download_shared_resources \
    --outputDir path/to/ref_dir
```

If you've already downloaded the shared resources, replace `--download_shared_resources` with `--shared_resources path/to/shared_resources` to avoid re-downloading these large files.

If `<genome name>` (the value passed to `--genome`) contains `hg19`, `hg38`, `mm10`, or `mm39` (or an alias, e.g. `GRCh38`, `GRCm39`), RENEE_NF automatically picks up the matching Arriba fusion-calling reference files from the directory given by `--arriba_db_dir`. For any other genome, pass `--fusion_blacklist`, `--fusion_cytoband`, `--fusion_protdomain`, and `--fusion_known_fusions` explicitly.

### Running

Once the genome resources have been built, run the pipeline against your samples. Replace `path/to/ref_dir` with the output directory used to build the genome resources above, and `<genome name>` with the genome name chosen there.

```
renee_nf run \
    -profile <singularity | docker> \
    --mode <local | slurm> \
    --input input.csv \
    -c path/to/ref_dir/genome/<genome name>.config \
    --genome <genome name> \
    --shared_resources path/to/ref_dir/shared_resources \
    --fastq_screen_conf path/to/ref_dir/shared_resources/fastq_screen_db/fastq_screen_p1.conf \
    --fastq_screen_conf2 path/to/ref_dir/shared_resources/fastq_screen_db/fastq_screen_p2.conf \
    --outputDir path/to/output
```

`--shared_resources` alone covers the Kraken2 database lookup, but FastQ Screen's two config files have no such fallback -- without `--fastq_screen_conf`/`--fastq_screen_conf2` set, FastQ Screen is silently skipped.

## Help & Contributing

Come across a **bug**? Open an [issue](https://github.com/CCBR/RENEE_NF/issues) and include a minimal reproducible example.

Have a **question**? Ask it in [discussions](https://github.com/CCBR/RENEE_NF/discussions).

Want to **contribute** to this project? Check out the [contributing guidelines](docs/contributing.md).

## References

This repo was originally generated from the [CCBR Nextflow Template](https://github.com/CCBR/CCBR_NextflowTemplate).
The template takes inspiration from nektool[^1] and the nf-core template.
If you plan to contribute your pipeline to nf-core, don't use this template -- instead follow nf-core's instructions[^2].

[^1]: nektool https://github.com/beardymcjohnface/nektool

[^2]: instructions for nf-core pipelines https://nf-co.re/docs/contributing/tutorials/creating_with_nf_core

[^3]: See also our reusable modules and subworkflows for CCBR nextflow pipelines: <https://github.com/CCBR/nf-modules>
