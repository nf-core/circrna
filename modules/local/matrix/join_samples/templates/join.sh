#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PY'

import polars as pl
from pathlib import Path

samples = "${samples.join(' ')}".split(" ")
matrices = "${matrices}".split(" ")
metacols = "${metacols}".split(",")

matrix_by_sample = {
    sample: next(matrix for matrix in matrices if Path(matrix).name.startswith(f"{sample}."))
    for sample in samples
}

dfs = {
    sample:
    pl.scan_csv(matrix,
                separator="\\t",
                new_columns=metacols + [sample],
                has_header=True)
        .with_columns(pl.col(sample).cast(pl.Float64, strict=False))
        .group_by(metacols).agg(pl.sum(sample).alias(sample))
    for sample, matrix in sorted(matrix_by_sample.items())}

df_order = pl.concat([df.select(metacols) for df in dfs.values()]).unique().sort(metacols)

dfs_sorted = [
    df_order.join(df, on=metacols, how="left", coalesce=True)
        .sort(metacols)
        .select(sample)
    for sample, df in dfs.items()
    ]

df = pl.concat([df_order] + dfs_sorted, how="horizontal").fill_null(0)

df.collect().write_csv("${meta.id}.joined.tsv", separator="\\t")

PY
