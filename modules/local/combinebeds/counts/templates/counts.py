#!/usr/bin/env python

import platform
import base64
import json

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
consider_strand = "${consider_strand}" == "true"
aggregation = "${aggregation}"
meta_id = "${meta.id}"
prefix = "${prefix}"
suffix = "${suffix}"

candidate_path = "${candidates}"
bed_paths = "${beds}".split()

columns = ["chr", "start", "end", "name", "score", "strand"]

df_candidates = pl.scan_csv(candidate_path, has_header=False, separator="\\t", new_columns=columns)
df_candidates = df_candidates.select(columns)
df_candidates = df_candidates.with_columns(sample=pl.lit("candidate"), tool=pl.lit("candidate"), score=pl.lit(None))

df = pl.scan_csv(bed_paths, has_header=False, separator="\\t", new_columns=columns + ["sample", "tool"])
df_combined = pl.concat([df, df_candidates])

df_combined = df_combined.sort("end"  ).with_columns(end_group  =pl.col("end"  ).diff().fill_null(0).gt(max_shift).cum_sum())
df_combined = df_combined.sort("start").with_columns(start_group=pl.col("start").diff().fill_null(0).gt(max_shift).cum_sum())

df_candidates = df_combined.filter(pl.col("sample") == "candidate")
df = df_combined.filter(pl.col("sample") != "candidate")

group_cols = ["chr", "start_group", "end_group"] + (["strand"] if consider_strand else [])
df = df.join(df_candidates, on=group_cols, how="inner")
df = df.filter((pl.col("start") - pl.col("start_right")).abs() <= max_shift)
df = df.filter((pl.col("end") - pl.col("end_right")).abs() <= max_shift)

df = df.select(["chr", "start", "end", "strand", "start_group", "end_group", "sample", "tool", "score"])

df = df.group_by(["chr", "strand", "start_group", "end_group", "start", "end"]).len().join(df, on=group_cols, how="inner")

df = df.filter((pl.col("start") - pl.col("start_right")).abs() <= max_shift)
df = df.filter((pl.col("end") - pl.col("end_right")).abs() <= max_shift)
df = df.group_by(["chr", "start", "end", "strand", "start_group", "end_group", "sample", "tool"]).agg(score=pl.sum("score"))
df = df.collect().lazy()

samples = df.select("sample").group_by("sample").len().collect()["sample"].to_list()
df = df.collect().pivot(on="sample", values="score", index=["chr", "start", "end", "strand", "start_group", "end_group"], aggregate_function=aggregation).lazy()
df = df.group_by(["chr", "strand", "start_group", "end_group"] + samples).agg(start=pl.col("start").first(), end=pl.col("end").first())
df = df.with_columns(id=pl.col("chr") + pl.lit(":") + pl.col("start").cast(str) + pl.lit("-") + pl.col("end").cast(str) + pl.lit(":") + pl.col("strand"))
df = df.sort("chr", "start", "end", "strand")
df = df.select(["id"] + samples)
df = df.fill_null(0)

df.sink_csv(f"{prefix}.{suffix}", separator="\\t", include_header=True)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
