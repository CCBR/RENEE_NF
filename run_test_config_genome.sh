#! /bin/bash
# Run the test config file

# load these modules if needed
# module load nextflow
# module load singularity
# samplesheet_invalidchars.csv
nextflow run -profile biowulf,test,slurm main.nf \
 -c conf/local_singularity.config \
 --input assets/samplesheet_single_read.csv \
 --genome hg38_36 \
 -resume
