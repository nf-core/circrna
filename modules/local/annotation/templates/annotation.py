#!/usr/bin/env python

import pandas as pd
import numpy as np
import platform
import csv

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


columns = {
    0: 'chr',
    1: 'start',
    2: 'end',
    3: 'name',
    4: 'score',
    5: 'strand',
    9: 'tx_start',
    10: 'tx_end',
    14: 'attributes'
}

attributes = ['gene_id', 'gene_name', 'transcript_id']

exon_boundary = int("${exon_boundary}")

df = pd.read_csv("${gtf_intersection}", sep="\\t", header=None, usecols=columns.keys())
df = df.rename(columns=columns)

# Extract circRNAs without match
mask = df['tx_start'] == -1
df_intergenic = df[mask]
df = df[~mask]
df_intergenic['type'] = 'intergenic-circRNA'
df_intergenic['gene_id'] = 'intergenic_' + df_intergenic['name']
df_intergenic['gene_name'] = 'intergenic_' + df_intergenic['name']
df_intergenic['transcript_id'] = 'intergenic_' + df_intergenic['name']

# Convert attributes to a dictionary
df['attributes'] = df['attributes'].apply(lambda row: dict([[value.strip(r'"') for value in entry.strip().split(' ', 1)] for entry in row.split(';') if entry]))
# Keep only the attributes we want
df['attributes'] = df['attributes'].apply(lambda row: {key: row[key] for key in attributes if key in row})
# Convert attributes to columns
df = pd.concat([df.drop(['attributes'], axis=1), df['attributes'].apply(pd.Series)], axis=1)

df['any_outside'] = (df['start'] < df['tx_start'] - exon_boundary) | (df['end'] > df['tx_end'] + exon_boundary)
# Perfect is inverse of any_outside
df['perfect'] = ~df['any_outside']
# Drop any_outside
df = df.drop(['any_outside', 'tx_start', 'tx_end'], axis=1)

df = df.groupby(['chr', 'start', 'end', 'strand']).aggregate({
    'name': lambda x: x.iloc[0],
    'score': lambda x: x.iloc[0],
    'gene_id': lambda x: list(x),
    'gene_name': lambda x: list(x),
    'transcript_id': lambda x: list(x),
    'perfect': lambda x: list(x)
})

def filter_perfect(row, col):
    if any(row['perfect']):
        matching_values = [value for value, perfectness in zip(row[col], row['perfect']) if perfectness]
    else:
        matching_values = row[col]
    valid_values = set([value for value in matching_values if type(value) == str])
    return ",".join(valid_values) if valid_values else "NaN"

def determine_type(row):
    if row["no_transcript"]:
        return "ciRNA"
    if any(row['perfect']):
        return "circRNA"
    else:
        return 'EI-circRNA'

df['no_transcript'] = df['transcript_id'].apply(lambda x: all([type(value) != str and np.isnan(value) for value in x]))
df['type'] = df.apply(lambda row: determine_type(row), axis=1)
df['gene_id'] = df.apply(lambda row: filter_perfect(row, 'gene_id'), axis=1)
df['gene_name'] = df.apply(lambda row: filter_perfect(row, 'gene_name'), axis=1)
df['transcript_id'] = df.apply(lambda row: filter_perfect(row, 'transcript_id'), axis=1)
# Drop perfect
df = df.drop(['perfect'], axis=1)

df = df.reset_index()
df_intergenic = df_intergenic.reset_index()
bed_order = ['chr', 'start', 'end', 'name', 'score', 'strand', 'type', 'gene_id', 'gene_name', 'transcript_id']
df = df[bed_order]
df_intergenic = df_intergenic[bed_order]

df = pd.concat([df, df_intergenic], axis=0)

db_intersections = "${db_intersections}".split()
has_db = len(db_intersections) > 0

if has_db:
    db_colnames = ['chr', 'start', 'end', 'name', 'score', 'strand', 'db_chr', 'db_start', 'db_end', 'db_name', 'db_score', 'db_strand']
    db_usecols = ['chr', 'start', 'end', 'name', 'score', 'strand', 'db_name']
    df_databases = pd.concat([pd.read_csv(db_path, sep="\\t", names=db_colnames, usecols=db_usecols) for db_path in db_intersections])

    # Group by chr, start, end, name, score, strand, and aggregate the db_name to list
    df_databases = df_databases.groupby(['chr', 'start', 'end', 'name', 'score', 'strand']).aggregate({
        'db_name': lambda x: ",".join([val for val in x if val != '.'])
    })

    df_databases['db_name'] = df_databases['db_name'].apply(lambda x: x if x else '.')

    df = df.merge(df_databases, how='left', on=['chr', 'start', 'end', 'name', 'score', 'strand'])
else:
    df['db_name'] = "."

# Sort by chr, start, end
df = df.sort_values(['chr', 'start', 'end'])

df.to_csv("${prefix}.bed", sep='\\t', index=False, header=False)

# Convert to GTF
df['source'] = 'circRNA'
df['frame'] = '.'
df['attributes'] = 'gene_id "' + df['gene_id'] + '"; gene_name "' + df['gene_name'] + '"; transcript_id "circ_' + df['name'] + '"; db_ids "' + df['db_name'] + '";'

gtf_order = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
df = df[gtf_order]

df.to_csv("${prefix}.gtf", sep='\\t', index=False, header=False, quoting=csv.QUOTE_NONE)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
