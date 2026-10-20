#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PYTHON'


import polars as pl

# Parameters

max_shift = int("${max_shift}")

# Logic

df_consensus = pl.read_csv("${consensus}", separator="\\t", has_header=False)
df_consensus = df_consensus.select("column_0", "column_1", "column_2")
df_consensus = df_consensus.rename({"column_0": "Chr", "column_1": "Start", "column_2": "End"})

df_consensus = df_consensus.with_columns(
    Start=pl.col("Start").cast(pl.Int64) + 1,
    End=pl.col("End").cast(pl.Int64),
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
df_anno = df_anno.with_columns(
    pl.col("Start").cast(pl.Int64),
    pl.col("End").cast(pl.Int64),
)
df_anno = df_anno.unique()

df_joined = df_consensus.join(df_anno, on=["Chr", "Start", "End"], how="inner")
df_joined.select("BSJ").collect().write_csv("${prefix}.txt", include_header=False)
PYTHON
