#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PYTHON'
# Fill unstranded BED records from canonical genomic BSJ splice motifs.

import pysam


def infer_strand(reference, fields):
    # Return a strand only for canonical flanks around a valid BED interval.
    if len(fields) < 6 or fields[5] != ".":
        return fields

    try:
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if start < 2 or end < start:
            return fields
        acceptor = reference.fetch(chrom, start - 2, start).upper()
        donor = reference.fetch(chrom, end, end + 2).upper()
    except (KeyError, ValueError, OSError):
        return fields

    if acceptor == "AG" and donor == "GT":
        fields[5] = "+"
    elif acceptor == "AC" and donor == "CT":
        fields[5] = "-"
    return fields


with pysam.FastaFile("${fasta}") as reference, open("${bed}") as source, open("${prefix}.bed", "w") as destination:
    for line in source:
        fields = line.rstrip("\\n").split("\\t")
        destination.write("\\t".join(infer_strand(reference, fields)) + "\\n")
PYTHON
