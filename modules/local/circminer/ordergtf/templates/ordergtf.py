#!/usr/bin/env python3

import platform
import sys
from itertools import groupby

import gffutils
from gffutils.iterators import DataIterator

exon_count = 0
transcript = None

with open("${prefix}.gtf", "w") as output:
    # AGAT groups exons by transcript in genomic order; buffer only that group.
    for feature_type, features in groupby(DataIterator("${gtf}"), key=lambda feature: feature.featuretype):
        if feature_type == "exon":
            if transcript is None:
                sys.exit("circMiner GTF ordering: exon feature occurred before a transcript feature")
            exons = list(features)
            exon_count += len(exons)
            output.writelines(str(exon) + "\\n" for exon in (reversed(exons) if transcript.strand == "-" else exons))
        elif feature_type in ("gene", "transcript"):
            for feature in features:
                output.write(str(feature) + "\\n")
                transcript = feature if feature_type == "transcript" else None

if not exon_count:
    sys.exit("circMiner GTF ordering: no exon features were found")
with open("versions.yml", "w", encoding="utf-8") as handle:
    print("${task.process}:", file=handle)
    print(f'    python: "{platform.python_version()}"', file=handle)
    print(f'    gffutils: "{gffutils.__version__}"', file=handle)
