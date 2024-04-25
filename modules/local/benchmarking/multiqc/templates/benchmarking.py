#!/usr/bin/env python3

import pandas as pd
import json
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

df = pd.read_csv("$jaccard", sep='\\t', index_col=0)

data = {
    "id": "benchmarking",
    "section_name": "Benchmarking",
    "description": "Benchmarking of the tools",
    "plot_type": "bargraph",
    "data": df.to_dict()
}

with open("benchmarking_mqc.json", "w") as f:
    json.dump(data, f)

versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}
with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))