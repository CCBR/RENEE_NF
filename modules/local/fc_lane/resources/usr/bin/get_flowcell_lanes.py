#!/usr/bin/env python3

"""Report read count, flowcell/lane identifiers, and MD5 for a FASTQ file."""

import gzip
import hashlib
import sys


def usage(message="", exitcode=0):
    """Display command usage and exit."""
    print(
        "Usage: python3 {} sample.R1.fastq.gz sample_name".format(sys.argv[0]),
        file=sys.stderr,
    )
    if message:
        print(message, file=sys.stderr)
    raise SystemExit(exitcode)


def get_flowcell_lane(sequence_identifier):
    """Return flowcell and lane values from common FASTQ identifier formats."""
    identifier_fields = sequence_identifier.strip().split(":")

    if len(identifier_fields) >= 7:
        # CASAVA >= 1.8:
        # @J00170:88:HNYVJBBXX:8:1101:6390:1244 1:N:0:ACTTGA
        return identifier_fields[2], identifier_fields[3]

    if sequence_identifier.startswith("@SRR"):
        # SRA identifiers may retain instrument/lane metadata:
        # @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
        try:
            return identifier_fields[0].split()[1], identifier_fields[1]
        except IndexError:
            # Otherwise use the accession as both the best available flowcell
            # and lane identifiers:
            # @SRR6755966.1 1 length=101
            accession = identifier_fields[0].split()[0].split(".")[0]
            return accession, accession.lstrip("@")

    # CASAVA < 1.8:
    # @HWUSI-EAS100R:6:73:941:1973#0/1
    try:
        return identifier_fields[0], identifier_fields[1]
    except IndexError as error:
        raise ValueError(
            "Unable to determine flowcell/lane from FASTQ identifier: "
            + sequence_identifier.strip()
        ) from error


def md5sum(filename, blocksize=65536):
    """Return the MD5 checksum of a file without loading it into memory."""
    hasher = hashlib.md5()
    with open(filename, "rb") as file_handle:
        for block in iter(lambda: file_handle.read(blocksize), b""):
            hasher.update(block)
    return hasher.hexdigest()


def main(filename, sample):
    """Write a two-line FASTQ information table to standard output."""
    flowcells = set()
    lanes = set()
    flowcell_lanes = set()
    line_count = 0

    with gzip.open(filename, "rt") as file_handle:
        for line_count, line in enumerate(file_handle, start=1):
            if (line_count - 1) % 4 == 0:
                flowcell, lane = get_flowcell_lane(line)
                flowcell = flowcell.lstrip("@")
                flowcells.add(flowcell)
                lanes.add(lane)
                flowcell_lanes.add("{}_{}".format(flowcell, lane))

    if line_count % 4:
        raise ValueError(
            "FASTQ contains {} lines; expected a multiple of four".format(line_count)
        )

    print(
        "sample_name\ttotal_read_pairs\tflowcell_ids\tlanes"
        "\tflowcell_lanes\tmd5_checksum"
    )
    print(
        "{}\t{}\t{}\t{}\t{}\t{}".format(
            sample,
            line_count // 4,
            ",".join(sorted(flowcells)),
            ",".join(sorted(lanes)),
            ",".join(sorted(flowcell_lanes)),
            md5sum(filename),
        )
    )


if __name__ == "__main__":
    if any(option in sys.argv for option in ("-h", "--help", "-help")):
        usage()
    if len(sys.argv) != 3:
        usage(
            message="Error: failed to provide all required positional arguments!",
            exitcode=1,
        )
    main(sys.argv[1], sys.argv[2])
