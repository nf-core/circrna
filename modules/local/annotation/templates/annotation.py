#!/usr/bin/env python

import platform
import csv

import pandas as pd
import numpy as np

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
    8: 'feature_type',
    9: 'feature_start',
    10: 'feature_end',
    14: 'attributes'
}

attributes = ['gene_id', 'gene_name', 'transcript_id']
feature_types = ['gene', 'exon']

try:
    df = pd.read_csv("${gtf_intersection}", sep="\\t", header=None, usecols=columns.keys())
except pd.errors.EmptyDataError:
    raise ValueError("Intersection between circRNAs and GTF file is empty.")
df = df.rename(columns=columns)
df 

# Extract circRNAs without match
mask = df['feature_start'] == -1
df_intergenic = df[mask]
df = df[~mask]
df_intergenic['type'] = 'intergenic-circRNA'
df_intergenic['gene_id'] = 'intergenic_' + df_intergenic['name']
df_intergenic['gene_name'] = 'intergenic_' + df_intergenic['name']
df_intergenic['transcript_id'] = 'intergenic_' + df_intergenic['name']

# Convert attributes to a dictionary
df['attributes'] = df['attributes'].apply(lambda row: dict([[value.strip(r'"') for value in entry.strip().split('=', 1)] for entry in row.split(';') if entry]))
# Make sure all attributes are present
df_incomplete = df['attributes'].apply(lambda row: ", ".join([key for key in attributes if key not in row]))
df_incomplete = df_incomplete[df_incomplete != ""]
if len(df_incomplete) > 0:
    counts = df_incomplete.value_counts()
    counts.name = 'count'
    counts.index.name = 'missing'
    raise ValueError(f"The following attributes are missing in the intersection file:\\n\\n{counts.to_frame()}")
# Keep only the attributes we want
df['attributes'] = df['attributes'].apply(lambda row: {key: row[key] for key in attributes if key in row})
# Convert attributes to columns
df = pd.concat([df.drop(['attributes'], axis=1), df['attributes'].apply(pd.Series)], axis=1)

df['contained'] = (df['start'] >= df['feature_start']) & (df['end'] <= df['feature_end'])
df = df.drop(['feature_start', 'feature_end'], axis=1)

df = df.groupby(['chr', 'start', 'end', 'strand']).aggregate({
    'name': 'first',
    'score': 'first',
    'gene_id': list,
    'gene_name': list,
    'transcript_id': list,
    'feature_type': list,
    'contained': list
})

def determine_type(row):
    if not any(row['contained']):
        return "partially_intergenic-circRNA"

    if "exon" in row["feature_type"]:
        if "intron" in row["feature_type"]:
            return "EI-circRNA"
        else:
            return "circRNA"
    if "intron" in row["feature_type"]:
        return "ciRNA"
    
    return "unknown-circRNA"

def get_representation(row, column):
    values = set(row[column])
    if len(values) == 0:
        return "."
    return ",".join(values)

df['type'] = df.apply(lambda row: determine_type(row), axis=1)
# Use df['name'] if gene or transcript column is empty
df['gene_id'] = df.apply(lambda row: get_representation(row, 'gene_id'), axis=1)
df['gene_name'] = df.apply(lambda row: get_representation(row, 'gene_name'), axis=1)
df['transcript_id'] = df.apply(lambda row: get_representation(row, 'transcript_id'), axis=1)
# Drop contained
df = df.drop(['contained'], axis=1)

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
