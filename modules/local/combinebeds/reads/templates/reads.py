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
meta_id = "${meta.id}"
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

bsj_columns = ["chr", "start", "end"] + (["strand"] if consider_strand else [])

# Investigate read-level agreement
df_bsj_groups = df.group_by(*bsj_columns, "reads").agg(pl.col("tool").unique()).collect()

# Plot read-level agreement
bsj_upset_data = from_memberships(df_bsj_groups["tool"].to_list())
bsj_upset = UpSet(bsj_upset_data, show_counts=True, subset_size="count")
bsj_upset.plot()
plot_file = "${prefix}:read_agreement.png"
plt.savefig(plot_file)
plt.close()

# MultiQC
image_string = base64.b64encode(open(plot_file, "rb").read()).decode("utf-8")
image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

multiqc = {
    'id': f"{meta_id}_read_agreement",
    'parent_id': "read_agreement",
    'parent_name': 'Read agreement',
    'parent_description': 'While different tools may agree on the presence of a BSJ, they may disagree on the reads that support it. These plots show the read-level agreement between tools for each BSJ.',
    'section_name': meta_id,
    'description': f'UpSet plot showing the read-level agreement between tools for each BSJ in {meta_id}',
    'plot_type': 'image',
    'data': image_html
}

with open(f"{meta_id}_read_agreement_mqc.json", "w") as f:
    f.write(json.dumps(multiqc, indent=4))

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
