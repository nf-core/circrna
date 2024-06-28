#!/usr/bin/env python3
import os
import re
import numpy as np

meta = "$meta"
depth = "$depth"
bed = "$bed"

def calculate_correlation(bed_file_path, depth_file_path):
    # Read the BED file and store circRNA regions
    circRNA_regions = []
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            parts = line.strip().split()
            chr = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            circRNA_regions.append((chr, start, end))

    # Read the sequence depth file
    circRNA_presence = []
    sequencing_depth = []

    with open(depth_file_path, 'r') as depth_file:
        bed_idx = 0
        bed_len = len(circRNA_regions)

        for line in depth_file:
            parts = line.strip().split()
            chr = parts[0]
            pos = int(parts[1])
            depth = int(parts[2])

            # Move bed_idx to the correct region
            while bed_idx < bed_len and (circRNA_regions[bed_idx][0] < chr or (circRNA_regions[bed_idx][0] == chr and circRNA_regions[bed_idx][2] < pos)):
                bed_idx += 1

            # Check if the current position is within a circRNA region
            if bed_idx < bed_len and circRNA_regions[bed_idx][0] == chr and circRNA_regions[bed_idx][1] <= pos <= circRNA_regions[bed_idx][2]:
                circRNA_presence.append(1)
            else:
                circRNA_presence.append(0)

            sequencing_depth.append(depth)

    # Calculate the correlation coefficient
    correlation = np.corrcoef(sequencing_depth, circRNA_presence)[0, 1]

    return correlation

id = re.search(r'id:(\\w+)', meta).group(1)

with open(depth, 'r') as paths:
    for path in paths:
        path = path.strip()
        basename = os.path.splitext(os.path.basename(path))[0]

        if basename == id:
            # Calculate correlation using the provided bed and depth files
            corr = calculate_correlation(bed, path)
            # Write the correlation to stats.txt
            with open('corr_mqc.tsv', 'w') as outfile:
                header = "tool\\tpearson_corr\\n"
                outfile.write(header + meta.split(",")[4][:-1] + '\\t' + str(corr))
