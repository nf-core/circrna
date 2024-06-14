#!/usr/bin/env python3

import platform
import pandas as pd

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

duplicates_fun = "${duplicates_fun}"

if duplicates_fun not in ['sum', 'mean', 'max', 'min']:
    raise ValueError(f"Invalid value for duplicates_fun: {duplicates_fun}")

columns = ['chr', 'start', 'end', 'name', 'count', 'strand']
dfs = [pd.read_csv(bed, sep='\t', header=None, names=columns) for bed in "${beds}".split()]
df = pd.concat(dfs)

df['tool_count'] = 1
df = df.groupby(['chr', 'start', 'end', 'strand', 'name']).agg({'count': duplicates_fun,
                                                        'tool_count': 'sum'}).reset_index()
df = df[df['tool_count'] >= int("${tool_filter}")]

df = df[columns]
df.to_csv("${meta.id}.merged.bed", sep='\t', index=False, header=False)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
