#!/usr/bin/env python3

"""
Parse the log file from preseq

Usage:
    python3 parse_preseq_log.py <path/to/preseq.log>

adapted from https://github.com/CCBR/Dockers/blob/4f245f4e418dd1908b7d91d487e45470d07096df/chipseq/ccbr_preseq/ccbr_nrf.py#L1-L17
"""


def main(log_filename):
    with open(log_filename, "r") as infile:
        for line in infile:
            if line.startswith("TOTAL READS"):
                tot_reads = float(line.strip().split("= ")[1])
            elif line.startswith("DISTINCT READS"):
                distinct_reads = float(line.strip().split("= ")[1])
            elif line.startswith("1\t"):
                one_pair = float(line.strip().split()[1])
            elif line.startswith("2\t"):
                two_pair = float(line.strip().split()[1])
    nrf = round(distinct_reads / tot_reads, 3)
    pbc1 = round(one_pair / distinct_reads, 3)
    pbc2 = round(one_pair / two_pair, 3)
    print(f"{nrf}\t{pbc1}\t{pbc2}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1])
