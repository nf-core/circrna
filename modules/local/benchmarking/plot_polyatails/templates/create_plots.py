#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import platform

# Function to calculate the average nA value for a given file
def calculate_average_na(file_path):
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t')  # Use tab as the separator
    
    # Calculate the mean of the nA column
    avg_na = df['nA'].mean()
    
    # Extract the sample name from the 'sample' column (since all rows have the same sample)
    sample_name = df['sample'].iloc[0]
    
    return avg_na, sample_name

def main(file_list, output_file):
    # Initialize lists to store the averages and corresponding samples
    averages = []
    samples = []
    
    # Loop through each file in the list
    for file_path in file_list:
        avg_na, sample_name = calculate_average_na(file_path)
        averages.append(avg_na)
        samples.append(sample_name)
    
    # Write the results to a TSV file
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')  # Tab-separated output
        writer.writerow(['sample', 'avg'])  # Write header
        for sample, avg in zip(samples, averages):
            writer.writerow([sample, avg])

if __name__ == "__main__":
    file_list_r = '$real'.split(' ')
    output_file_r = 'real.tsv'
    main(file_list_r, output_file_r)
    file_list_b = '$benchmarking'.split(' ')
    output_file_b = 'benchmarking.tsv'
    main(file_list_b, output_file_b)

    