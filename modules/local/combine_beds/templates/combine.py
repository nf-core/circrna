#!/usr/bin/env python

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

max_shift = int("${max_shift}")
min_files = int("${min_files}")

df = pl.scan_csv("*.bed", 
                 separator="\\t",
                 has_header=False,
                 include_file_paths="file",
                 new_columns=["chr", "start", "end", "name", "score", "strand"])

for col in ["end", "start"]:
    df = df.sort(["chr", col])
    df = df.with_columns(**{f"{col}_group": pl.col(col).diff().fill_null(0).gt(max_shift).cum_sum()})

df = (df.group_by(["chr", "start_group", "end_group"])
    .agg(pl.col("start").median().round().cast(int),
         pl.col("end").median().round().cast(int),
         pl.col("file").n_unique().alias("n_files"))
    .filter(pl.col("n_files") >= min_files)
    .select(["chr", "start", "end"])
    .with_columns(name=pl.col("chr").cast(str) + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str),
                  score=pl.lit("."),
                  strand=pl.lit(".")))

df.collect().write_csv("${prefix}.${suffix}", separator="\\t", include_header=False)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
