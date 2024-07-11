#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

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
