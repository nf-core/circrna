#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Merge counts across tools for a single sample')
parser.add_argument('--beds', type=str, nargs='+', help='List of bed files to merge')
parser.add_argument('--tool_filter', type=int, help='Minimum number of tools to keep a circRNA')
parser.add_argument('--duplicates_fun', type=str, help='Function to apply to duplicates', choices=['sum', 'mean', 'max', 'min'])
parser.add_argument('--output', type=str, help='Output file')

args = parser.parse_args()

columns = ['chr', 'start', 'end', 'strand', 'count']
dfs = [pd.read_csv(bed, sep='\t', header=None, names=columns) for bed in args.beds]
df = pd.concat(dfs)

df['tool_count'] = 1

df = df.groupby(['chr', 'start', 'end', 'strand']).agg({'count': args.duplicates_fun,
                                                        'tool_count': 'sum'}).reset_index()
df = df[df['tool_count'] >= args.tool_filter]

df.drop('tool_count', axis=1, inplace=True)
df["name"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str) + ":" + df["strand"]

df = df[['chr', 'start', 'end', 'name', 'count', 'strand']]

df.to_csv(args.output, sep='\t', index=False, header=False)
