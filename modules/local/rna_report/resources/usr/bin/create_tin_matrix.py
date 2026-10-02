#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Combines per-sample RSeQC TIN (*.tin.xls) files into a single transcript-by-sample
matrix. Ported from the RENEE Snakemake pipeline's create_tin_matrix.py so both
pipelines resolve duplicate transcript IDs (e.g. PAR genes annotated on both
chrX and chrY) the same way: last value wins, one row per transcript ID.
"""

from __future__ import print_function
import sys
import os
import pandas


def create(file, tin_dict, key_index=0, parse_index=4):
    """Populates the TIN nested dictionary
    @param file <str>: Path to RSEQC output file with TIN values to extract
    @param tin_dict <dict>: Dictionary to populate where [samplebasename][transcriptid] = tin_value
    @param key_index <int>: Index of the field to join multiple files
    @param parse_index <int>: Index of field of interest (i.e. TIN value)
    """

    with open(file, "r") as fh:
        header = next(fh).strip().split("\t")
        colid = header[key_index]
        sample = os.path.basename(file).split(".tin.xls")[0]

        for line in fh:
            linelist = line.strip().split("\t")
            tid = linelist[key_index]
            tinvalue = linelist[parse_index]
            if sample not in tin_dict:
                tin_dict[sample] = {}

            tin_dict[sample][tid] = tinvalue

    return colid, tin_dict


if __name__ == "__main__":
    # Get filenames to parse
    args = sys.argv
    files = sys.argv[1:]

    # Check if at least one file was provided
    if not len(args) >= 2:
        print("FATAL: Failed to provide at least one input file!")
        sys.exit("Usage:\n python {} *.tin.xls > combined_TIN.tsv".format(args[0]))

    # Populate tins with TIN values for all transcripts across all samples
    tins = {}
    for file in files:
        keycolname, tins = create(file, tins)

    df = pandas.DataFrame(tins)
    # Print dataframe to standard output
    df.to_csv(sys.stdout, sep="\t", header=True, index=True, index_label=keycolname)
