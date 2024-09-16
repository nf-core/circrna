#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
import platform

input_bed_file_1 = '$bedfile1'
input_bed_file_2 = '$bedfile2'

# Read input files
def read_bed_file(file_path, label):
    data = {'chromosome': [], 'start': [], 'strand': [], 'file_label': []}
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\\t')
        for row in reader:
            data['chromosome'].append(row[0])
            data['start'].append(int(row[1]))
            data['strand'].append(row[3])
            data['file_label'].append(label.lower())
    return data

def combine_data(data1, data2):
    for key in data1:
        data1[key].extend(data2[key])
    return data1

data1 = read_bed_file(input_bed_file_1, 'total')
data2 = read_bed_file(input_bed_file_2, 'polya')

# Combine the two datasets
combined_data = combine_data(data1, data2)

# Create a DataFrame
df = pd.DataFrame({
    'Chromosome': combined_data['chromosome'],
    'Start Location': combined_data['start'],
    'Strand': combined_data['strand'],
    'File Label': combined_data['file_label']
})

# Sort DataFrame to ensure consistent plotting order
df.sort_values(by=['File Label', 'Strand'], inplace=True)

# Plotting
fig, ax = plt.subplots(figsize=(12, 6))

palette = {
    "total +": "red",
    "total -": "lightcoral",
    "polya +": "blue",
    "polya -": "lightblue"
}

# Draw violins
for file_label in df['File Label'].unique():
    sns.violinplot(
        x="Chromosome",
        y="Start Location",
        hue="Strand",
        data=df[df['File Label'] == file_label],
        palette={"+" : palette[f"{file_label} +"], "-" : palette[f"{file_label} -"]},
        split=True,
        ax=ax,
        scale="count",
        scale_hue=False,
        saturation=0.75,
        inner=None
    )

# Set transparency for all violins
for violin in ax.collections:
    violin.set_alpha(0.25)

# Legend
custom_lines = [
    Line2D([0], [0], color=palette[f"{file} {strand}"], lw=4, alpha=0.25)
    for file in df['File Label'].unique()
    for strand in ["+", "-"]
]
ax.legend(
    custom_lines,
    [f"{file} : {strand}" for file in df['File Label'].unique() for strand in ["+", "-"]],
    title="File : Strand"
)

plt.title('Start Locations of circRNA by Chromosome and Strand')

plot_file_name = f"{input_bed_file_1.replace('.bed','')}_{input_bed_file_2.replace('.bed','')}_mqc.png"

# Save the plot
plt.savefig(plot_file_name, bbox_inches='tight')

#version capture
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

versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "seaborn": sns.__version__
    }
}
with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
