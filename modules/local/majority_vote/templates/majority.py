#!/usr/bin/env python3 
import pandas as pd
import platform


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


df = pd.read_csv("$bindingsites",
                          sep="\\t", header=0, names=['mirna', 'target', 'start', 'end', 'tool' ])
df = df.groupby(['mirna', 'target'])['tool'].apply(set).reset_index()

# performing majority vote keeping only mirna binding sites that meet the required number of votes
min_tools = int("${params.mirna_vote}")
df = df[df['tool'].apply(len) >= min_tools].copy()
df = df.drop('tool', axis=1)

df.to_csv('${meta.id}.majority.tsv', sep='\\t', index=False)

df = df.groupby('target')['mirna'].apply(lambda x: ','.join(x)).reset_index()
df.to_csv('${meta.id}.targets.tsv', sep='\\t', index=False, header=False)

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
