#!/usr/bin/env python3

import pandas as pd
import platform

tsv = "$tsv"
data = pd.read_csv(tsv, sep='\\t')

data['tool'] = data['tool'].str.replace('tool:', '')
average_values = data.groupby('tool')['pearson_corr'].mean().reset_index()

output_file_path = tsv
average_values.to_csv(output_file_path, sep='\\t', index=False)

#version capture
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

versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}
with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))

