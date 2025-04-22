#!/usr/bin/env python

import platform
import base64
import json
import yaml
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import upsetplot
from upsetplot import from_memberships, UpSet

import polars as pl

tools = "${tools.join(" ")}".split()
bed_files = "${beds}".split()
consider_strand = "${consider_strand}" == "true"
min_tools = int("${min_tools}")
dfs = []

for path, tool in zip(bed_files, tools):
    df = pl.scan_csv(path, separator="\\t", has_header=False)
    df = df.with_columns(tool=pl.lit(tool))
    dfs.append(df)

df = pl.concat(dfs, rechunk=True)
columns = ["chr", "start", "end", "name", "score", "strand", "reads"]
df = df.rename({f"column_{i+1}": col for i, col in enumerate(columns)})
df = df.with_columns(
    reads = pl.col("reads").str.split(",")
)

df = df.explode("reads")

# Investigate read-level agreement
df_read_groups = df.group_by("reads").agg(pl.col("tool").unique()).collect()

# Plot read-level agreement
read_upset_data = from_memberships(df_read_groups["tool"].to_list())
read_upset = UpSet(read_upset_data, show_counts=True, subset_size="count")
read_upset.plot()
plt.savefig("${prefix}:read_agreement.png")
plt.close()

bsj_columns = ["chr", "start", "end"] + (["strand"] if consider_strand else [])

# Investigate BSJ-level agreement
df_bsj_groups = df.group_by(*bsj_columns, "reads").agg(pl.col("tool").unique()).collect()

# Plot BSJ-level agreement
bsj_upset_data = from_memberships(df_bsj_groups["tool"].to_list())
bsj_upset = UpSet(bsj_upset_data, show_counts=True, subset_size="count")
bsj_upset.plot()
plt.savefig("${prefix}:bsj_agreement.png")
plt.close()

df_bsj_groups = df_bsj_groups.with_columns(
    tool_count = pl.col("tool").list.len()
)

df_filtered = df_bsj_groups.filter(pl.col("tool_count") >= min_tools)
df_filtered = df_filtered.group_by(bsj_columns).agg(pl.col("reads").unique())
df_filtered = df_filtered.with_columns(
    reads = pl.col("reads").list.join(","),
    score = pl.col("reads").list.len()
)

# Name should look like FUSIONJUNC_27/13 where 27 is the row number and 13 is the score
df_filtered = df_filtered.with_columns(
    name = pl.lit("FUSIONJUNC_") + pl.arange(0, df_filtered.height).cast(str) + pl.lit("/") + pl.col("score").cast(str)
)

if not consider_strand:
    df_filtered = df_filtered.with_columns(
        strand = pl.lit(".")
    )

df_filtered = df_filtered.select("chr", "start", "end", "name", "score", "strand", "reads")
df_filtered.write_csv("${prefix}.bed", separator="\\t", include_header=False)

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
        "seaborn": sns.__version__,
        "matplotlib": matplotlib.__version__,
        "upsetplot": upsetplot.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
