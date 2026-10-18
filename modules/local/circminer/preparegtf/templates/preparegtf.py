#!/usr/bin/env python3

import gzip
import platform
import re
import sys
from collections import OrderedDict


ATTRIBUTE = re.compile(r'\s*([^\s;]+)\s+"((?:[^"\\]|\\.)*)"\s*')


def fail(message):
    raise ValueError("circMiner GTF preparation: " + message)


def parse_attributes(text, line_number):
    attributes = OrderedDict()
    position = 0
    while position < len(text):
        match = ATTRIBUTE.match(text, position)
        if not match:
            fail("line {} has malformed GTF attributes".format(line_number))
        key, value = match.groups()
        if key in attributes:
            fail("line {} repeats attribute '{}'".format(line_number, key))
        attributes[key] = value
        position = match.end()
        if position == len(text) or text[position] != ';':
            fail("line {} has malformed GTF attributes".format(line_number))
        position += 1
    return attributes


def format_attributes(attributes):
    return " ".join('{} "{}";'.format(key, value) for key, value in attributes.items())


def open_gtf(filename):
    with open(filename, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    return gzip.open(filename, "rt") if compressed else open(filename)


genes = OrderedDict()
transcript_to_gene = {}
exon_count = 0

try:
    with open_gtf("${gtf}") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                fail("line {} does not have 9 tab-separated fields".format(line_number))
            feature = fields[2]
            if feature != "exon":
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                fail("line {} has non-integer coordinates".format(line_number))
            if start < 1 or end < start:
                fail("line {} has invalid coordinates {}-{}".format(line_number, start, end))
            if fields[6] not in ("+", "-"):
                fail("line {} has invalid strand '{}'".format(line_number, fields[6]))
            attributes = parse_attributes(fields[8], line_number)
            gene_id = attributes.get("gene_id")
            transcript_id = attributes.get("transcript_id")
            if not gene_id or not transcript_id:
                fail("line {} exon requires non-empty gene_id and transcript_id attributes".format(line_number))

            gene = genes.setdefault(gene_id, {"chromosome": fields[0], "strand": fields[6], "transcripts": OrderedDict(), "name": attributes.get("gene_name")})
            if gene["chromosome"] != fields[0] or gene["strand"] != fields[6]:
                fail("gene_id '{}' is used on multiple chromosomes or strands".format(gene_id))
            if attributes.get("gene_name") and not gene["name"]:
                gene["name"] = attributes["gene_name"]
            owner = transcript_to_gene.setdefault(transcript_id, gene_id)
            if owner != gene_id:
                fail("transcript_id '{}' is used by multiple genes".format(transcript_id))
            transcript = gene["transcripts"].setdefault(transcript_id, {"exons": [], "name": attributes.get("gene_name")})
            transcript["exons"].append((start, end, fields, attributes))
            exon_count += 1
except (OSError, ValueError) as error:
    sys.exit(str(error))

if not exon_count:
    sys.exit("circMiner GTF preparation: no exon features were found")

with open("${prefix}.gtf", "w") as output:
    ordered_genes = sorted(genes.items(), key=lambda item: (item[1]["chromosome"], min(exon[0] for transcript in item[1]["transcripts"].values() for exon in transcript["exons"]), item[0]))
    for gene_id, gene in ordered_genes:
        transcripts = gene["transcripts"]
        gene_exons = [exon for transcript in transcripts.values() for exon in transcript["exons"]]
        gene_attributes = OrderedDict([("gene_id", gene_id)])
        if gene["name"]:
            gene_attributes["gene_name"] = gene["name"]
        output.write("\t".join([gene["chromosome"], "circminer_preparegtf", "gene", str(min(exon[0] for exon in gene_exons)), str(max(exon[1] for exon in gene_exons)), ".", gene["strand"], ".", format_attributes(gene_attributes)]) + "\n")
        ordered_transcripts = sorted(transcripts.items(), key=lambda item: (min(exon[0] for exon in item[1]["exons"]), item[0]))
        for transcript_id, transcript in ordered_transcripts:
            exons = sorted(transcript["exons"], key=lambda exon: (exon[0], exon[1]), reverse=gene["strand"] == "-")
            transcript_attributes = OrderedDict([("gene_id", gene_id), ("transcript_id", transcript_id)])
            if transcript["name"] or gene["name"]:
                transcript_attributes["gene_name"] = transcript["name"] or gene["name"]
            output.write("\t".join([gene["chromosome"], "circminer_preparegtf", "transcript", str(min(exon[0] for exon in exons)), str(max(exon[1] for exon in exons)), ".", gene["strand"], ".", format_attributes(transcript_attributes)]) + "\n")
            for exon_number, (start, end, fields, attributes) in enumerate(exons, 1):
                attributes = OrderedDict((key, value) for key, value in attributes.items() if key != "exon_number")
                attributes["exon_number"] = str(exon_number)
                fields[8] = format_attributes(attributes)
                output.write("\t".join(fields) + "\n")

with open("versions.yml", "w") as versions:
    versions.write('"${task.process}":\n    python: "{}"\n'.format(platform.python_version()))
