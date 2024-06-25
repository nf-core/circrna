#!/usr/bin/env python3

import pandas as pd

tsv = "$tsv"
data = pd.read_csv(tsv, sep='\t')

data['tool'] = data['tool'].str.replace('tool:', '')
average_values = data.groupby('tool')['pearson_corr'].mean().reset_index()

output_file_path = tsv
average_values.to_csv(output_file_path, sep='\t', index=False)

