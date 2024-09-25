#!/usr/bin/env python3

import platform

import polars as pl
import yaml

samples = "${samples.join(' ')}".split(" ")
matrices = "${matrices}".split(" ")
metacols = "${metacols}".split(",")

dfs = {
    sample:
    pl.scan_csv(matrix,
                separator="\\t",
                new_columns=metacols + [sample],
                has_header="${has_header}" == "true")
        .group_by(metacols).agg(pl.sum(sample).alias(sample))
    for sample, matrix in zip(samples, matrices)}

df_order = pl.concat([df.select(metacols) for df in dfs.values()]).unique().sort(metacols[0])

dfs_sorted = [
    df_order.join(df, on=metacols, how="left", coalesce=True)
        .select(sample)
    for sample, df in dfs.items()
    ]

df = pl.concat([df_order] + dfs_sorted, how="horizontal").fill_null(0)

df.collect().write_csv("${meta.id}.joined.tsv", separator="\\t")

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "polars": pl.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))
