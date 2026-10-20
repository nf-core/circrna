#!/usr/bin/env python3

import platform
import importlib.metadata

# Versions
with open("versions.yml", "w") as handle:
    handle.write('"${task.process}":\\n')
    handle.write(f'  python: "{platform.python_version()}"\\n')
    handle.write(f'  polars: "{importlib.metadata.version("polars")}"\\n')

import polars as pl

df_annotation = pl.scan_csv(
    "${bsj_annotation}",
    separator="\\t",
    has_header=False,
    new_columns=["chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "exonCount", "exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron"]
)
df_reads = pl.scan_csv(
    "${bsj_reads}",
    separator="\\t",
    has_header=False,
    new_columns=["chr", "start", "end", "name", "score", "strand", "reads"]
)

join_cols = ["chr", "start", "end", "strand"]

df_merged = df_reads.join(
    df_annotation.select(join_cols + ["geneName"]),
    on=join_cols,
    how="inner"
)
df_merged = df_merged.select(
    "chr", "start", "end", "geneName", "strand", "reads"
).collect()

df_merged.write_csv("${prefix}.txt", separator="\\t", include_header=False)
