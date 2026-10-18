#!/usr/bin/env python3
"""Convert psiRC candidate calls and resolve strand from supporting transcripts."""

import re

import yaml


def transcript_strands(transcriptome):
    strands = {}
    with open(transcriptome) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            # Upstream psiRC emits: chr:start-end_transcript_Exon_Lengths_..._Offsets_..._+/-
            match = re.match(r"^>.*?:\d+-\d+_(.+)_Exon_Lengths_[^_]+_Offsets_[^_]+_([+-])\s*$", line)
            if match:
                strands[match.group(1)] = match.group(2)
    return strands


def supporting_strand(name, strands):
    support = set()
    for junction in name.split("="):
        match = re.match(r"(.+)_E\d+B\d+$", junction)
        if match and match.group(1) in strands:
            support.add(strands[match.group(1)])
    return support.pop() if len(support) == 1 else "."


strands = transcript_strands("${transcriptome}")
with open("${bed}") as source, open("${prefix}.bed", "w") as destination:
    for number, line in enumerate(source):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        try:
            score, names = fields[3].split(":", 1)
            int(score)
            int(fields[1])
            int(fields[2])
        except ValueError:
            continue
        destination.write("\t".join([fields[0], fields[1], fields[2], f"FUSIONJUNC_{number}/{score}", score, supporting_strand(names.split(":", 1)[0], strands)]) + "\n")

with open("versions.yml", "w") as handle:
    yaml.safe_dump({"${task.process}": {"python": "3"}}, handle)
