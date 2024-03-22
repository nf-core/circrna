#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Prepare matrices for circTest')
parser.add_argument('--in_genes', type=str, help='Gene input tsv file')
parser.add_argument('--in_circs', type=str, help='CircRNA input tsv file')
parser.add_argument('--out_genes', type=str, help='Gene output tsv file')
parser.add_argument('--out_circs', type=str, help='CircRNA output tsv file')

args = parser.parse_args()

df_genes = pd.read_csv(args.in_genes, sep='\t', index_col=0)
df_circs = pd.read_csv(args.in_circs, sep='\t', index_col=0)

gene_ids = df_circs["gene_id"]
df_circs.to_csv(args.out_circs, sep='\t')

df_genes = df_genes.loc[gene_ids]
df_genes.to_csv(args.out_genes, sep='\t')