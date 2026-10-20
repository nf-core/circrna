#!/usr/bin/env python3

import platform
import importlib.metadata

# Versions
with open("versions.yml", "w") as handle:
    handle.write('"${task.process}":\\n')
    handle.write(f'  python: "{platform.python_version()}"\\n')
    handle.write(f'  polars: "{importlib.metadata.version("polars")}"\\n')

import polars as pl

paths = "${bindingsites}".split(" ")

df = pl.scan_csv(paths,
                 separator="\\t",
                 has_header=False,
                 new_columns=['mirna', 'target', 'start', 'end', 'tool'])

df = df.select(["mirna", "target", "tool"])

df = df.group_by(['mirna', 'target']).agg(pl.col("tool").n_unique())

df = df.filter(pl.col("tool") >= int("${min_tools}")) \
    .select(["mirna", "target"])

df = df.collect()

df.write_csv('${meta.id}.majority.tsv', separator='\\t', include_header=False)

# Create targets file

df = df.group_by('mirna').agg(pl.col("target").str.join(","))

df.write_csv('${meta.id}.targets.tsv', separator='\\t', include_header=False)
