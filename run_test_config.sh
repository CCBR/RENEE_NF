#! /bin/bash
# Run the test config file

# load these modules if needed
# module load nextflow
# module load singularity
nextflow run -profile test,biowulf,slurm main.nf
