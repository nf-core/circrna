#!/usr/bin/env python3

import platform
import sys


exon_count = 0
transcript_exons = []
transcript_strand = None

with open("${gtf}") as source, open("${prefix}.gtf", "w") as output:
    for line in source:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\\t")
        feature = fields[2]

        if feature == "gene":
            if transcript_exons:
                output.writelines(reversed(transcript_exons) if transcript_strand == "-" else transcript_exons)
                transcript_exons = []
            output.write(line)
            transcript_strand = None
        elif feature == "transcript":
            if transcript_exons:
                output.writelines(reversed(transcript_exons) if transcript_strand == "-" else transcript_exons)
                transcript_exons = []
            output.write(line)
            transcript_strand = fields[6]
        elif feature == "exon":
            if transcript_strand is None:
                sys.exit("circMiner GTF ordering: exon feature occurred before a transcript feature")
            transcript_exons.append(line)
            exon_count += 1

    if transcript_exons:
        output.writelines(reversed(transcript_exons) if transcript_strand == "-" else transcript_exons)

if not exon_count:
    sys.exit("circMiner GTF ordering: no exon features were found")

with open("versions.yml", "w") as versions:
    versions.write('"${task.process}":\\n    python: "{}"\\n'.format(platform.python_version()))
