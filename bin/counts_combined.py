#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Merge counts across tools for a single sample')
parser.add_argument('--beds', type=str, nargs='+', help='List of bed files to merge')
parser.add_argument('--out_bed', type=str, help='Output bed file')
parser.add_argument('--out_tsv', type=str, help='Output tsv file')

args = parser.parse_args()

columns = ['chr', 'start', 'end', 'name', 'count', 'strand']
dfs = {os.path.basename(bed).split('.')[0]: pd.read_csv(bed,
                   sep='\t',
                   header=None,
                   index_col=["chr", "start", "end", "strand"],
                   usecols=["chr", "start", "end", "strand", "count"],
                   names=columns) for bed in args.beds}

dfs = [df.rename(columns={'count': sample}) for sample, df in dfs.items()]
df = pd.concat(dfs, axis=1)
df = df.fillna(0)
df.to_csv(args.out_bed, sep='\t', header=True, index=True)

df.index = df.index.map(lambda x: f'{x[0]}:{x[1]}-{x[2]}:{x[3]}')
df.index.name = 'ID'
df.to_csv(args.out_tsv, sep='\t')
