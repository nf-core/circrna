#!/usr/bin/env python

import platform

import polars as pl
import yaml

# Parameters

max_shift = int("${max_shift}")

# Logic

df_consensus = pl.read_csv("${consensus}", separator="\\t", has_header=False)
df_consensus = df_consensus.select("column_1", "column_2", "column_3")
df_consensus = df_consensus.rename({"column_1": "Chr", "column_2": "Start", "column_3": "End"})

df_consensus = df_consensus.with_columns(
    Start=pl.col("Start") + 1
)

df_length = df_consensus.shape[0]
df_consensus = df_consensus.lazy()

shifts = list(range(-max_shift, max_shift + 1))
df_consensus = df_consensus.with_columns(
    start_shifts=pl.repeat(shifts, df_length),
    end_shifts=pl.repeat(shifts, df_length)
)

df_consensus = df_consensus.explode("start_shifts").explode("end_shifts")
df_consensus = df_consensus.with_columns(
    Start=pl.col("Start") + pl.col("start_shifts"),
    End=pl.col("End") + pl.col("end_shifts")
)
df_consensus = df_consensus.select("Chr", "Start", "End")
df_consensus = df_consensus.unique()

df_anno = pl.scan_csv("${ciri_annotation}".split(), separator="\\t", has_header=True, truncate_ragged_lines=True)
df_anno = df_anno.select("BSJ", "Chr", "Start", "End")
df_anno = df_anno.unique()

df_joined = df_consensus.join(df_anno, on=["Chr", "Start", "End"], how="inner")
df_joined.select("BSJ").collect().write_csv("${prefix}.txt", include_header=False)

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
