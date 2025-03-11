#! /usr/bin/env python3

import platform

import polars as pl

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

df_counts = pl.scan_csv("${counts}", separator="\\t")
df_counts = df_counts.rename({"paired.junctions": "Count"})

join_columns = ['Chr', 'Start', 'End']
df_coordinates = pl.scan_csv("${coordinates}", separator="\\t")
df_coordinates = df_coordinates.select(join_columns + ['Strand'])

df_counts = df_counts.join(df_coordinates, on=join_columns, how='left')

df_counts = df_counts.collect()
df_counts = df_counts.with_columns(
    name=pl.lit("FUSIONJUNC_") + pl.int_range(0, len(df_counts)).cast(pl.Utf8) + pl.lit("/") + pl.col("Count").cast(pl.Utf8)
)
df_counts = df_counts.select(["Chr", "Start", "End", "name", "Count", "Strand"])
df_counts.write_csv("${prefix}.${suffix}", separator="\\t")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
