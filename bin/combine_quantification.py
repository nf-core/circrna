#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Merge quantification files into a single file')
parser.add_argument('--inputs', type=str, nargs='+', help='List input files to merge')
parser.add_argument('--output', type=str, help='Path to output combined quantification file')

args = parser.parse_args()

dfs = [pd.read_csv(f, sep='\t', index_col=[0, 1]) for f in args.inputs]
df_combined = pd.concat(dfs, axis=1, sort=True)

df_combined.to_csv(args.output, sep='\t')
