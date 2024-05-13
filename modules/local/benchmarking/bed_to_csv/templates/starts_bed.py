import pandas as pd
import sys

# Command line arguments
input_bed_file_1 = sys.argv[1]
output_csv_file_1 = sys.argv[2]
input_bed_file_2 = sys.argv[3]
output_csv_file_2 = sys.argv[4]

# Process the first BED file
data1 = pd.read_csv(input_bed_file_1, sep='\t', header=None, usecols=[0, 1, 5], names=['chromosome', 'start', 'strand'])
data1.to_csv(output_csv_file_1, index=False)

# Process the second BED file
data2 = pd.read_csv(input_bed_file_2, sep='\t', header=None, usecols=[0, 1, 5], names=['chromosome', 'start', 'strand'])
data2.to_csv(output_csv_file_2, index=False)