#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PY'

import base64
import json

import polars as pl
import upsetplot
import matplotlib
import matplotlib.pyplot as plt


# Parameters

max_shift = int("${max_shift}")
consider_strand = "${consider_strand}" == "true"
min_tools = int("${min_tools}")
min_samples = int("${min_samples}")
meta_id = "${meta.id}"
prefix = "${prefix}"

# Logic

df = pl.scan_csv("*.bed",
                 separator="\\t",
                 has_header=False,
                 new_columns=["chr", "start", "end", "name", "score", "strand", "sample", "tool"])

df = df.group_by("chr", "start", "end", "strand").agg(tools=pl.col("tool").unique(), samples=pl.col("sample").unique())

df = df.sort("end"  ).with_columns(end_group  =pl.col("end"  ).diff().fill_null(0).gt(max_shift).cum_sum())
df = df.sort("start").with_columns(start_group=pl.col("start").diff().fill_null(0).gt(max_shift).cum_sum())

group_cols = ["chr", "start_group", "end_group"] + (["strand"] if consider_strand else [])
df = df.join(df, on=group_cols, how="inner")
df = df.filter((pl.col("start") - pl.col("start_right")).abs() <= max_shift)
df = df.filter((pl.col("end") - pl.col("end_right")).abs() <= max_shift)

df = df.select(["chr", "start", "end", "strand", "samples_right", "tools_right"])
df = df.group_by(["chr", "start", "end", "strand"]).agg(**{
    "samples": pl.col("samples_right").list.explode().unique(),
    "tools": pl.col("tools_right").list.explode().unique()
}).with_columns(n_samples=pl.col("samples").map_elements(lambda x: len(x), return_dtype=int),
                n_tools=pl.col("tools").map_elements(lambda x: len(x), return_dtype=int))

df = df.collect()
df = df.with_columns(name=pl.lit("FUSIONJUNC_") + pl.arange(0, len(df)).cast(str) + pl.lit("/0"),
                    score=pl.lit("."))

df_filtered = df.filter(
    (pl.col("n_tools") >= min_tools) & (pl.col("n_samples") >= min_samples)
).select("chr", "start", "end", "name", "score", "strand")

if df_filtered.is_empty():
    exit(0)

df_filtered.write_csv("${prefix}.${suffix}", separator="\\t", include_header=False)

# Plots

for col in ["samples", "tools"]:
    memberships = df[col].to_list()
    if len({member for group in memberships for member in group}) <= 1:
        continue
    dataset = upsetplot.from_memberships(memberships)

    n_combinations = 50
    if len(dataset) > n_combinations:
        dataset = dataset.nlargest(n_combinations)

    upsetplot.plot(dataset,
                   orientation="horizontal",
                   show_counts=True,
                   subset_size="count")
    # UpSetPlot uses one-element arrays for count-label margins; Matplotlib
    # requires scalar coordinates with current NumPy releases.
    for axis in plt.gcf().axes:
        for label in axis.texts:
            label.set_position(tuple(value.item() if hasattr(value, "item") else value
                                     for value in label.get_position()))
    plot_file = f"{prefix}_{col}.upset.png"
    plt.savefig(plot_file)
    plt.close()

    with open(plot_file, "rb") as image:
        image_string = base64.b64encode(image.read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'
    multiqc = {
        "id": f"{meta_id}_upset_{col}",
        "parent_id": "upset_plots",
        "parent_name": "UpSet Plots",
        "parent_description": "UpSet plots showing the overlap between tools for each sample",
        "section_name": f"UpSet {col}: {meta_id} ",
        "description": f"UpSet plot showing the overlap between {col} for {meta_id}",
        "plot_type": "image",
        "data": image_html,
    }
    with open(f"{prefix}_{col}.upset_mqc.json", "w") as handle:
        json.dump(multiqc, handle, indent=4)

PY
