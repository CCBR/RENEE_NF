#!/usr/bin/env python3
import sys

unknown = '"Unknown";'

# One row per transcript_id. Not every GTF has an explicit `transcript`
# feature line for every transcript -- e.g. NCBI/GenBank-derived GTFs only
# give `gene`/CDS/exon lines, and even a well-formed Ensembl GTF can have an
# orphan exon line with no matching `transcript` line if it's been trimmed
# down to a subset (as in subworkflows/local/prepare_genome/tests' fixture) --
# so derive the isoform entry from any line that has a transcript_id,
# preferring an explicit `transcript` line's values over a CDS/exon-derived
# guess when both exist for the same id.
isoforms = {}
order = []

for i in list(
    filter(
        lambda x: not x[0].startswith("#"),
        list(map(lambda x: x.strip().split("\t"), open(sys.argv[1]).readlines())),
    )
):
    gene_id = ""
    j = i[8].split()
    transcript_id = unknown
    gene_id = unknown
    gene_name = unknown
    transcript_name = unknown
    for k in list(range(0, len(j) - 1, 2)):
        if j[k] == "transcript_id":
            transcript_id = j[k + 1]
        elif j[k] == "gene_id":
            gene_id = j[k + 1]
        elif j[k] == "transcript_name":
            transcript_name = j[k + 1]
        elif j[k] == "gene_name":
            gene_name = j[k + 1]
    if transcript_id == unknown:
        continue
    if transcript_name == unknown and transcript_id != unknown:
        transcript_name = transcript_id
    if gene_name == unknown and gene_id != unknown:
        gene_name = gene_id

    is_transcript_line = i[2] == "transcript"
    seen = transcript_id in isoforms
    if not seen:
        order.append(transcript_id)
    if not seen or (is_transcript_line and not isoforms[transcript_id][3]):
        isoforms[transcript_id] = (
            transcript_name,
            gene_id,
            gene_name,
            is_transcript_line,
        )

for transcript_id in order:
    transcript_name, gene_id, gene_name, _ = isoforms[transcript_id]
    s = "%s  %s  %s  %s" % (
        transcript_id[:-1],
        transcript_name[:-1],
        gene_id[:-1],
        gene_name[:-1],
    )
    print(s)
