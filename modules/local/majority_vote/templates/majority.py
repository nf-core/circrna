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


predictions = pd.read_csv("$bindingsites",
                          sep="\\t", header=0, names=['mirna', 'target', 'start', 'end', 'tool' ])

# start = False  
# end = False 
# TODO: Ask if this means an exact match or an "either ..., or ..."
complete = True  
majority = ${params.mirna_vote}

# if start:  # group by start indices
#     predictions = predictions.groupby(['mirna', 'target', 'start'])['tool'].apply(set).reset_index()
# elif end:  # group by end indices
#     predictions = predictions.groupby(['mirna', 'target', 'end'])['tool'].apply(set).reset_index()
# elif complete:  # group by both indices
#     predictions = predictions.groupby(['mirna', 'target', 'start', 'end'])['tool'].apply(set).reset_index()

predictions = predictions.groupby(['mirna', 'target', 'start', 'end'])['tool'].apply(set).reset_index()

# performing majority vote keeping only mirna binding sites that meet the required number of votes
post_vote_predictions = predictions[predictions['tool'].apply(len) >= majority].copy()
out = post_vote_predictions.drop('tool', axis=1)


out.to_csv('${meta.id}.majority.tsv', sep='\\t', index=False)

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
