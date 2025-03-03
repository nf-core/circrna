#!/usr/bin/env python

import platform
import base64
import json

import polars as pl
import upsetplot
import matplotlib
import matplotlib.pyplot as plt

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

max_shift = int("${max_shift}")
consider_strand = "${consider_strand}" == "true"
min_tools = int("${min_tools}")
min_samples = int("${min_samples}")
meta_id = "${meta.id}"
prefix = "${prefix}"

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
    "samples": pl.col("samples_right").flatten().unique(),
    "tools": pl.col("tools_right").flatten().unique()
}).with_columns(n_samples=pl.col("samples").map_elements(lambda x: len(x), return_dtype=int),
                n_tools=pl.col("tools").map_elements(lambda x: len(x), return_dtype=int))

df = df.with_columns(name=pl.col("chr") + ":" + pl.col("start").cast(str) + "-" + pl.col("end").cast(str) + ":" + pl.col("strand"),
                    score=pl.lit("."))

df_aggregated = df.collect().to_pandas()
n_bsjs = len(df_aggregated)

df_filtered = df_aggregated[(df_aggregated["n_tools"] >= min_tools) & (df_aggregated["n_samples"] >= min_samples)]
df_filtered = df_filtered[["chr", "start", "end", "name", "score", "strand"]]

df_filtered.to_csv("${prefix}.${suffix}", sep="\\t", header=False, index=False)

for col in ["samples", "tools"]:
    series = df_aggregated[col]
    if series.explode().nunique() <= 1:
        continue
    memberships = series.to_list()
    dataset = upsetplot.from_memberships(memberships)
    upsetplot.plot(dataset,
                   orientation='horizontal',
                   show_counts=True,
                   subset_size="count",
                   min_degree=2,
                   min_subset_size=min(50, int(n_bsjs * 0.02)))
    plot_file = f"{prefix}_{col}.upset.png"
    plt.savefig(plot_file)

    image_string = base64.b64encode(open(plot_file, "rb").read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    multiqc = {
        'id': f"{meta_id}_upset_{col}",
        'parent_id': "upset_plots",
        'parent_name': 'UpSet Plots',
        'parent_description': 'UpSet plots showing the overlap between tools for each sample',
        'section_name': f'UpSet {col}: {meta_id} ',
        'description': f'UpSet plot showing the overlap between {col} for {meta_id}',
        'plot_type': 'image',
        'data': image_html
    }

    with open(f"{prefix}_{col}.upset_mqc.json", "w") as f:
        f.write(json.dumps(multiqc, indent=4))

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
        "upsetplot": upsetplot.__version__,
        "matplotlib": matplotlib.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
