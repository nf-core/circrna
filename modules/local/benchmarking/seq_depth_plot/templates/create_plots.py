#!/usr/bin/env python3
import logging
import numpy as np

# Initialize logging
logging.basicConfig(filename='process_log.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

meta = "$meta"
depth = "$depth"
bed = "$bed"

def extract_matching_paths(array, filepath):
    array_id = array.split(",")[0].split(":")[1]
    matching_paths = []
    with open(filepath, 'r') as file:
        for line in file:
            line_id = line.split(",")[0].split(":")[1]
            if array_id == line_id:
                path = line.split()[-1]
                matching_paths.append(path)
    return matching_paths

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

with open('stats.txt', 'w') as outfile:
    for path in extract_matching_paths(meta, depth):
        corr = calculate_correlation(bed, path)
        output_line = meta.split(",")[3] + '      ' + str(corr)
        outfile.write(output_line)
        logging.debug(f'Written to stats.txt: {output_line.strip()}')
