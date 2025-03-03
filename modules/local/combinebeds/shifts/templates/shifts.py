#!/usr/bin/env python

import platform
import base64
import json
from itertools import product

import polars as pl
import altair as alt

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

meta_id = "${meta.id}"

df = pl.scan_csv("${beds}".split(" "),
                    separator="\\t",
                    has_header=False,
                    new_columns=["chr", "start", "end", "name", "score", "strand", "sample", "tool"])

df = df.group_by("chr", "start", "end", "strand").agg(tools=pl.col("tool").unique(), samples=pl.col("sample").unique())

def get_group_sizes(df: pl.LazyFrame, max_shift: int, consider_strand: bool) -> pl.LazyFrame:
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

    df_sample_counts = (df.group_by("n_samples").agg(count=pl.col("n_samples").count())
        .sort("n_samples")
        .rename({"n_samples": "value"})
        .with_columns(metric=pl.lit("n_samples")))
    df_tool_counts = (df.group_by("n_tools").agg(count=pl.col("n_tools").count())
        .sort("n_tools")
        .rename({"n_tools": "value"})
        .with_columns(metric=pl.lit("n_tools")))

    df = (pl.concat([df_sample_counts, df_tool_counts])
        .with_columns(max_shift=max_shift, consider_strand=consider_strand))

    return df

shifts = [0, 1, 2, 3, 4, 5, 10, 20, 50]
consider_strands = [True, False]

dfs = []
for max_shift, consider_strand in product(shifts, consider_strands):
    df_ = get_group_sizes(df, max_shift, consider_strand)
    dfs.append(df_.collect())

df = pl.concat(dfs).with_columns(max_shift=pl.col("max_shift"))

metrics = {
    "n_samples": "Number of samples",
    "n_tools": "Number of tools"
}

for metric, title in metrics.items():
    n_unique = df.filter(pl.col("metric") == metric).select("value").n_unique("value")
    df_ = df.filter(pl.col("metric") == metric)
    plot = df_.plot.bar(x="max_shift:N", y="count", color=f"value:{("Q" if n_unique > 10 else "N")}", column="consider_strand")
    plot = plot.properties(title=title)

    plot_file = f"{metric}.png"
    plot.save(plot_file)

    image_string = base64.b64encode(open(plot_file, "rb").read()).decode("utf-8")
    image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

    multiqc = {
        'id': f"{meta_id}_shifts_{metric}",
        'parent_id': "shift_plots",
        'parent_name': 'Shift Plots',
        'parent_description': 'Stacked bar plots showing the agreement between tools and samples for different shift values',
        'section_name': title,
        'description': f'Stacked bar plot showing the agreement between tools and samples for different shift values, {"considering" if consider_strand else "ignoring"} strand',
        'plot_type': 'image',
        'data': image_html
    }

    with open(f"{metric}.shifts_mqc.json", "w") as f:
        f.write(json.dumps(multiqc, indent=4))

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
        "altair": alt.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
