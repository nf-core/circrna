#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd
import platform

def read_bed(bed_file: str, label: str):
    """Read BED file and return as DataFrame with an additional label column for color coding."""
    df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    df['label'] = label
    return df

def plot_bed(datasets: list, output_file: str):
    """Plot BED intervals on a genomic track with color coding.

    Args:
        datasets (list of tuples): List containing tuples of (DataFrame, color).
        output_file (str): Path to save the plot image.
    """
    fig, ax = plt.subplots(figsize=(10, 1))
    colors = {'real': 'blue', 'benchmarking': 'red'}

    for data, label in datasets:
        chromosome = data['chr'].unique()[0]
        intervals = [(row['start'], row['end'] - row['start']) for index, row in data.iterrows()]
        yrange = (0, 10)
        color = colors[label]
        ax.add_collection(BrokenBarHCollection(intervals, yrange, facecolors=color, label=label))

    ax.set_yticks([])
    ax.set_xlim(min(data['start'].min() for data, label in datasets), max(data['end'].max() for data, label in datasets))
    ax.set_title("Genomic loci by source")
    ax.legend(title="Source")
    plt.savefig(output_file)

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string."""
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\n"
    return yaml_str

def save_versions():
    """Save Python and pandas versions in YAML-like format."""
    versions = {
        "plot_bed.py": {
            "python": platform.python_version(),
            "pandas": pd.__version__
        }
    }
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))

real_bed_file = sys.argv[1]
benchmarking_bed_file = sys.argv[2]
output_file = sys.argv[3]

real_data = read_bed(real_bed_file, 'real')
benchmarking_data = read_bed(benchmarking_bed_file, 'benchmarking')

datasets = [(real_data, 'real'), (benchmarking_data, 'benchmarking')]
plot_bed(datasets, output_file)
save_versions()
