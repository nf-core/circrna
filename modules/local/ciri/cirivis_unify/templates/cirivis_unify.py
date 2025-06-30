#!/usr/bin/env python3
"""
Convert CIRI-vis output to BED12 format
"""

import pandas as pd
import yaml

def simple_bed12_conversion(df):
    """Simple conversion to BED12 format"""
    bed12_data = []
    
    for _, row in df.iterrows():
        # Parse exon coordinates (skip 0-0 breaks)
        exons = [coord for coord in row['isoform_cirexon'].split(',') 
                if coord.strip() and coord.strip() != '0-0']
        
        if not exons:
            continue
            
        # Calculate blocks
        block_sizes = []
        block_starts = []
        
        for exon in exons:
            start, end = map(int, exon.split('-'))
            block_sizes.append(end - start)
            block_starts.append(start - row['start'])  # Relative to chromStart
        
        # Create BED12 row
        bed12_data.append([
            row['Chr'],                                    # chrom
            row['start'],                                  # chromStart  
            row['end'],                                    # chromEnd
            row['Circle_ID'],                              # name
            row['isoform_exp'],                           # score (just use raw expression)
            row['strain'],                                 # strand
            row['start'],                                  # thickStart
            row['end'],                                    # thickEnd
            '0',                                          # itemRgb (black)
            len(block_sizes),                             # blockCount
            ','.join(map(str, block_sizes)) + ',',        # blockSizes
            ','.join(map(str, block_starts)) + ','        # blockStarts
        ])
    
    return pd.DataFrame(bed12_data, columns=[
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'
    ])

# Hardcoded paths
input_file = "${list}"
output_file = "${prefix}.bed"

print(f"Reading CIRI-vis data from: {input_file}")

# Read CIRI-vis output
df = pd.read_csv(input_file, sep='\\t')
print(f"Loaded {len(df)} circRNA entries")

# Convert to BED12
df_bed12 = simple_bed12_conversion(df)
print(f"Converted {len(df_bed12)} entries to BED12 format")

# Save BED12 file
df_bed12.to_csv(output_file, sep='\\t', header=False, index=False)
print(f"Saved BED12 file: {output_file}")

# Quick stats
print(f"\\nBlock count distribution:")
block_counts = df_bed12['blockCount'].value_counts().sort_index()
for blocks, count in block_counts.items():
    print(f"  {count} circRNAs with {blocks} block(s)")


# Versions

versions = {
    "${task.process}": {
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))
