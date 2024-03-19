#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Merge quantification files into a single file')
parser.add_argument('--inputs', type=str, nargs='+', help='List input files to merge')
parser.add_argument('--tx2gene', type=str, help='Path to tx2gene file')
parser.add_argument('--circ_annotation', type=str, help='Path to circRNA annotation file')
parser.add_argument('--out_linear', type=str, help='Path to output linear RNA quantification file')
parser.add_argument('--out_circular', type=str, help='Path to output circRNA quantification file')
parser.add_argument('--out_combined', type=str, help='Path to output combined quantification file')

args = parser.parse_args()

df_annotation = pd.read_csv(args.circ_annotation,
                            sep='\t', header=None,
                            names=['chr', 'start', 'end', 'name', 'score', 'strand', 'type', 'gene', 'transcript'],
                            index_col=['name'])

df_tx2gene = pd.read_csv(args.tx2gene, sep='\t')

dfs = [pd.read_csv(f, sep='\t', index_col=["gene_id", "gene_name"]) for f in args.inputs]
df_combined = pd.concat(dfs, axis=1, sort=True)

# Split linear and circular
linear_mask = df_combined.index.get_level_values(0) \
                .isin(df_tx2gene["gene_id"])
df_linear = df_combined[linear_mask].copy()
df_circular = df_combined[~linear_mask].copy()

df_linear.to_csv(args.out_linear, sep='\t')

# Perform annotation
df_annotation["annotation"] = df_annotation.apply(lambda row: 'intergenic' if row['type'] == 'intergenic-circRNA' else row['gene'], axis=1)
annotation_dict = df_annotation["annotation"].to_dict()

df_circular.reset_index(inplace=True)
df_circular["gene_name"] = df_circular['gene_id'].map(annotation_dict)
df_circular.set_index(["gene_id", "gene_name"], inplace=True)

df_circular.to_csv(args.out_circular, sep='\t')

# Save combined
df_linear["type"] = "gene"
df_circular["type"] = "circRNA"

df_combined = pd.concat([df_linear, df_circular], axis=0, sort=True)

df_combined.to_csv(args.out_combined, sep='\t')
