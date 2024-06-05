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

# Combine File Label and Strand into a single column for easier plotting
df['File_Strand'] = df['File Label'] + ' ' + df['Strand']

# Separate the data by strand
df_plus = df[df['Strand'] == '+']
df_minus = df[df['Strand'] == '-']

# Create figure and axes
fig, ax = plt.subplots(figsize=(12, 6))

# Define the color palette
palette_plus = {
    "real +": "red", 
    "benchmark +": "blue"
}
palette_minus = {
    "real -": "lightcoral",  
    "benchmark -": "lightblue"
}

# Draw violins for the + strand
sns.violinplot(
    x="Chromosome", 
    y="Start Location", 
    hue="File_Strand",
    data=df_plus,
    palette=palette_plus,
    split=True,
    ax=ax,
    scale="count",
    scale_hue=False,
    saturation=0.75,
    inner=None
)

# Draw violins for the - strand
sns.violinplot(
    x="Chromosome", 
    y="Start Location", 
    hue="File_Strand",
    data=df_minus,
    palette=palette_minus,
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
    Line2D([0], [0], color=color, lw=4, alpha=0.25) 
    for color in palette_plus.values()
] + [
    Line2D([0], [0], color=color, lw=4, alpha=0.25) 
    for color in palette_minus.values()
]
ax.legend(
    custom_lines, 
    list(palette_plus.keys()) + list(palette_minus.keys()), 
    title="File : Strand"
)

plt.title('Start Locations of circRNA by Chromosome and Strand')
plot_file_name = f"{input_bed_file_1.replace('.bed','')}_{input_bed_file_2.replace('.bed','')}_mqc.png"
# Save the plot
plt.savefig(plot_file_name, bbox_inches='tight')


