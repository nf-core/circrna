#!/usr/bin/env python3
import csv
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D

# Inputs directly from Nextflow variables
input_bed_file_1 = '$bedfile1'
input_bed_file_2 = '$bedfile2'

# Read input files function
def read_bed_file(file_path, label):
    data = {'chromosome': [], 'start': [], 'strand': [], 'file_label': []}
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            data['chromosome'].append(row[0])
            data['start'].append(int(row[1]))
            data['strand'].append(row[3])  # Assuming strand information is in the 4th column
            data['file_label'].append(label)  # Label to indicate the source file
    return data

# Combine data function
def combine_data(data1, data2):
    for key in data1:
        data1[key].extend(data2[key])
    return data1

data1 = read_bed_file(input_bed_file_1, 'real')
data2 = read_bed_file(input_bed_file_2, 'benchmark')

# Combine the two datasets
combined_data = combine_data(data1, data2)

# Create a DataFrame
df = pd.DataFrame({
    'Chromosome': combined_data['chromosome'],
    'Start Location': combined_data['start'],
    'Strand': combined_data['strand'],
    'File Label': combined_data['file_label']
})

# Separate the data by file label
df_real = df[df['File Label'] == 'real']
df_benchmark = df[df['File Label'] == 'benchmark']

# Create figure and axes
fig, ax = plt.subplots(figsize=(12, 6))

# Define the color palette
palette_real = {
    "+": "red",
    "-": "lightcoral"
}
palette_benchmark = {
    "+": "blue",
    "-": "lightblue"
}

# Draw violins for the real file
sns.violinplot(
    x="Chromosome",
    y="Start Location",
    hue="Strand",
    data=df_real,
    palette=palette_real,
    split=True,
    ax=ax,
    scale="count",
    scale_hue=False,
    saturation=0.75,
    inner=None
)

# Draw violins for the benchmark file
sns.violinplot(
    x="Chromosome",
    y="Start Location",
    hue="Strand",
    data=df_benchmark,
    palette=palette_benchmark,
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

# Compose a custom legend
custom_lines = [
    Line2D([0], [0], color=palette_real["+"], lw=4, alpha=0.25),
    Line2D([0], [0], color=palette_real["-"], lw=4, alpha=0.25),
    Line2D([0], [0], color=palette_benchmark["+"], lw=4, alpha=0.25),
    Line2D([0], [0], color=palette_benchmark["-"], lw=4, alpha=0.25)
]
ax.legend(custom_lines, ["Total +", "Total -", "PolyA +", "PolyA -"], title="File : Strand")

plt.title('Start Locations of circRNA by Chromosome and Strand')
plot_file_name = f"{input_bed_file_1.replace('.bed','')}_{input_bed_file_2.replace('.bed','')}_mqc.png"
# Save the plot
plt.savefig(plot_file_name, bbox_inches='tight')
