#!/usr/bin/env python

import platform
import yaml

import polars as pl

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))

# Main

max_shift = int("${max_shift}")
consider_strand = "${consider_strand}" == "true"
aggregation = "${aggregation}"
meta_id = "${meta.id}"
prefix = "${prefix}"
suffix = "${suffix}"

candidate_path = "${candidates}"
bed_paths = "${beds}".split()

columns = ["chr", "start", "end", "name", "score", "strand"]

try:
    df_candidates = pl.scan_csv(candidate_path, has_header=False, separator="\\t", new_columns=columns, raise_if_empty=True)
except pl.exceptions.NoDataError:
    print("No data in ${candidates}")
    with open(f"{prefix}.{suffix}", "w") as f:
        f.write('')
    exit(0)

df_candidates = df_candidates.select(columns)
df_candidates = df_candidates.with_columns(sample=pl.lit("candidate"), tool=pl.lit("candidate"), score=pl.lit(None))

try:
    df = pl.scan_csv(bed_paths, has_header=False, separator="\\t", new_columns=columns + ["sample", "tool"], raise_if_empty=True)
except pl.exceptions.NoDataError:
    print("No data in ${beds}")
    with open(f"{prefix}.{suffix}", "w") as f:
        f.write('')
    exit(0)

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

try:
    df = df.collect().lazy()
except pl.exceptions.NoDataError:
    print("No data after processing")
    with open(f"{prefix}.{suffix}", "w") as f:
        f.write('')
    exit(0)

samples = df.select("sample").group_by("sample").len().collect()["sample"].to_list()
df = df.collect().pivot(on="sample", values="score", index=["chr", "start", "end", "strand", "start_group", "end_group"], aggregate_function=aggregation).lazy()
df = df.group_by(["chr", "strand", "start_group", "end_group"] + samples).agg(start=pl.col("start").first(), end=pl.col("end").first())
df = df.with_columns(id=pl.col("chr") + pl.lit(":") + pl.col("start").cast(str) + pl.lit("-") + pl.col("end").cast(str) + pl.lit(":") + pl.col("strand"))
df = df.sort("chr", "start", "end", "strand")
df = df.select(["id"] + samples)
df = df.fill_null(0)

df.sink_csv(f"{prefix}.{suffix}", separator="\\t", include_header=True)
