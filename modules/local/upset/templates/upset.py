#!/usr/bin/env python3

import pandas as pd
import platform
import upsetplot
import matplotlib
import matplotlib.pyplot as plt
import distutils.version
import base64
import json

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

df_tools = pd.DataFrame(
    {
        "tool": "${tools.join(' ')}".split(" "),
        "file": "${beds.join(' ')}".split(" ")
    }
)

tool_files = df_tools.groupby("tool")["file"].apply(lambda x: set(x)).to_dict()
tool_ids = {}

for tool, files in tool_files.items():
    df_tool = pd.concat([pd.read_csv(f, sep="\\t", header=None) for f in files])
    tool_ids[tool] = set(df_tool[3].unique())

dataset = upsetplot.from_contents(tool_ids)

upsetplot.plot(dataset, orientation='horizontal', show_counts=True)
plot_file = "${meta.id}.upset.png"
plt.savefig(plot_file)

image_string = base64.b64encode(open(plot_file, "rb").read()).decode("utf-8")
image_html = f'<div class="mqc-custom-content-image"><img src="data:image/png;base64,{image_string}" /></div>'

multiqc = {
    'id': "${meta.id}_upset",
    'parent_id': "upset_plots",
    'parent_name': 'UpSet Plots',
    'parent_description': 'UpSet plots showing the overlap between tools for each sample',
    'section_name': 'UpSet: ${meta.id}',
    'description': 'UpSet plot showing the overlap between tools for sample ${meta.id}',
    'plot_type': 'image',
    'data': image_html
}

with open("${meta.id}.upset_mqc.json", "w") as f:
    f.write(json.dumps(multiqc, indent=4))

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "upsetplot": upsetplot.__version__,
        "matplotlib": matplotlib.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
