#!/usr/bin/env python3

import platform
import pandas as pd
import os

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

columns = ['chr', 'start', 'end', 'name', 'count', 'strand']
dfs = {os.path.basename(bed).split('.')[0]: pd.read_csv(bed,
                   sep='\t',
                   header=None,
                   index_col=["chr", "start", "end", "strand"],
                   usecols=["chr", "start", "end", "strand", "count"],
                   names=columns) for bed in "${beds}".split()}

dfs = [df.rename(columns={'count': sample}) for sample, df in dfs.items()]
df = pd.concat(dfs, axis=1)
df = df.fillna(0)
df.to_csv("merged_counts.bed", sep='\\t', header=True, index=True)

df.index = df.index.map(lambda x: f'{x[0]}:{x[1]}-{x[2]}:{x[3]}')
df.index.name = 'ID'
df.to_csv("merged_counts.tsv", sep='\\t')

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
