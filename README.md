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

### Build the genome and shared resources

Initialize the output directory, then run the build:

```
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

If `<genome name>` (the value passed to `--genome`) contains `hg19`, `hg38`, `mm10`, or `mm39` (or an alias, e.g. `GRCh38`, `GRCm39`), RENEE_NF automatically detects the matching Arriba fusion-calling reference files in the shared resources directory populated by `--download_shared_resources` above (or in the directory passed via `--arriba_db_dir`, if you set that instead). For any other genome, pass `--fusion_blacklist`, `--fusion_cytoband`, `--fusion_protdomain`, and `--fusion_known_fusions` explicitly.

### Building a second genome

When building a second genome, skip redownloading the large shared resources files by omitting `--download_shared_resources` and instead including `-c path/to/ref_dir/shared_resources.config`.

```
renee_nf run \
    -profile <singularity | docker> \
    -c path/to/ref_dir/shared_resources.config \
    --mode <local | slurm> \
    --build \
    --publish_dir_mode copy \
    --genome <second genome name> \
    --genome_fasta /path/to/second_genome.fa \
    --genes_gtf /path/to/second_genome.gtf
```

### Running

Once the genome resources have been built, run the pipeline against your samples. Replace `path/to/ref_dir` with the output directory used to build the genome resources above, and `<genome name>` with the genome name chosen there.

```
renee_nf run \
    -profile <singularity | docker> \
    --mode <local | slurm> \
    --input input.csv \
    --genome <genome name> \
    -c path/to/ref_dir/genome/<genome name>.config \
    -c path/to/ref_dir/shared_resources.config \
    --outputDir path/to/output
```

The `-c` flag lets you include additional config files, each supplying extra parameters for the pipeline. Outside of Biowulf, this pipeline needs two additional configs:

- `<genome name>.config` adds paths specific to the genome being used.
- `shared_resources.config` points to the Kraken2, FastQ Screen, and Arriba databases, which are shared across all genome builds.

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
