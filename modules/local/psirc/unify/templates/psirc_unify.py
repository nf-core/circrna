#!/usr/bin/env python3
"""
Convert PSIRC output to BED12 format
"""

import pandas as pd
import yaml

def psirc_to_bed12(df):
    """Convert PSIRC output to BED12 format"""
    bed12_data = []
    
    for _, row in df.iterrows():
        # Parse genomic_loci: "chrI:10046438-10047864"
        if pd.isna(row['genomic_loci']) or row['genomic_loci'] == '.':
            continue
            
        try:
            # Split chr and coordinates
            chr_part, coords = row['genomic_loci'].split(':')
            start, end = map(int, coords.split('-'))
        except:
            continue  # Skip malformed entries
        
        # Handle score - use back_splice_junction_read_count if available and numeric, otherwise skip row
        score = 0
        bsj = row['back_splice_junction_read_count']
        if pd.isna(bsj) or bsj == '.':
            continue
        try:
            score = int(float(bsj))
        except (ValueError, TypeError):
            continue  # Skip this row if not numeric

        # PSI-RC doesn't have detailed exon structure, so create single block
        block_size = end - start

        # Create BED12 row
        bed12_data.append([
            chr_part,                    # chrom
            start,                       # chromStart
            end,                         # chromEnd  
            row['output_name'],          # name
            score,                       # score
            row['strand'],               # strand
            start,                       # thickStart
            end,                         # thickEnd
            '0',                         # itemRgb
            1,                           # blockCount (single block)
            f"{block_size},",            # blockSizes
            "0,"                         # blockStarts (starts at 0)
        ])
    
    return pd.DataFrame(bed12_data, columns=[
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'
    ])

input_file = "${isoforms_tsv}"
output_file = "${prefix}.bed"

print(f"Reading PSIRC data from: {input_file}")

# Read PSIRC output
df = pd.read_csv(input_file, sep='\\t')
print(f"Loaded {len(df)} PSIRC entries")

# Convert to BED12
df_bed12 = psirc_to_bed12(df)
print(f"Converted {len(df_bed12)} entries to BED12 format")

# Save BED12 file
df_bed12.to_csv(output_file, sep='\\t', header=False, index=False)
print(f"Saved BED12 file: {output_file}")

# Quick stats
print(f"\\nPSIRC statistics:")
print(f"- Total entries: {len(df_bed12)}")
print(f"- Entries with read counts: {len(df_bed12[df_bed12['score'] > 0])}")
print(f"- Strand distribution:")
for strand, count in df_bed12['strand'].value_counts().items():
    print(f"  {strand}: {count}")

# Versions

versions = {
    "${task.process}": {
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))
