#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import platform

input_bed_file = '$bed'

# Read input file
def read_bed_file(file_path):
    data = {'chromosome': [], 'start': [], 'strand': []}
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            data['chromosome'].append(row[0])
            data['start'].append(int(row[1]))
            data['strand'].append(row[3])
    return data

# Read the data
data = read_bed_file(input_bed_file)

# Create a DataFrame
df = pd.DataFrame({
    'Chromosome': data['chromosome'],
    'Start Location': data['start'],
    'Strand': data['strand']
})

# Plotting
plt.figure(figsize=(12, 6))
sns.violinplot(x='Chromosome', y='Start Location', hue='Strand', data=df, split=True, inner="quartile")
plt.title('Start Locations of circRNA by Chromosome and Strand')
plt.legend(title="Strand")
plot_file_name = f"{input_bed_file.replace('.bed','')}_overlaps_mqc.png"

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
