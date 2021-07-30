#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    "nf-core/circrna": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "BBDUK": ["v_bbduk.txt", r"(\S+)"],
    "BEDTools": ["v_bedtools.txt", r"(\S+)"],
    "Bowtie": ["v_bowtie.txt", r"(\S+)"],
    "Bowtie2": ["v_bowtie2.txt", r"(\S+)"],
    "BWA": ["v_bwa.txt", r"(\S+)"],
    "CIRCexplorer2": ["v_circexplorer2.txt", r"(\S+)"],
    "CIRIquant": ["v_ciriquant.txt", r"(\S+)"],
    "Java": ["v_java.txt", r"(\S+)"],
    "MapSplice": ["v_mapsplice.txt", r"(\S+)"],
    "miRanda": ["v_miranda.txt", r"(\S+)"],
    "Perl": ["v_perl.txt", r"(\S+)"],
    "SamToFastq": ["v_picard.txt", r"(\S+)"],
    "Python": ["v_python.txt", r"(\S+)"],
    "R": ["v_R.txt", r"(\S+)"],
    "SAMtools": ["v_samtools.txt", r"(\S+)"],
    "Segemehl": ["v_segemehl.txt", r"(\S+)"],
    "STAR": ["v_star.txt", r"(\S+)"],
    "StringTie": ["v_stringtie.txt", r"(\S+)"],
}
results = OrderedDict()
results["nf-core/circrna"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["BBDUK"] = '<span style="color:#999999;">N/A</span>'
results["BEDtools"] = '<span style="color:#999999;">N/A</span>'
results["Bowtie"] = '<span style="color:#999999;">N/A</span>'
results["Bowtie2"] = '<span style="color:#999999;">N/A</span>'
results["BWA"] = '<span style="color:#999999;">N/A</span>'
results["CIRCexplorer2"] = '<span style="color:#999999;">N/A</span>'
results["CIRIquant"] = '<span style="color:#999999;">N/A</span>'
results["Java"] = '<span style="color:#999999;">N/A</span>'
results["MapSplice"] = '<span style="color:#999999;">N/A</span>'
results["miRanda"] = '<span style="color:#999999;">N/A</span>'
results["Perl"] = '<span style="color:#999999;">N/A</span>'
results["SamToFastq"] = '<span style="color:#999999;">N/A</span>'
results["Python"] = '<span style="color:#999999;">N/A</span>'
results["R"] = '<span style="color:#999999;">N/A</span>'
results["SAMtools"] = '<span style="color:#999999;">N/A</span>'
results["STAR"] = '<span style="color:#999999;">N/A</span>'
results["StringTie"] = '<span style="color:#999999;">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/circrna Software Versions'
section_href: 'https://github.com/nf-core/circrna'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
