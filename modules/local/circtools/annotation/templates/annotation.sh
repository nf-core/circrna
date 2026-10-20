#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PYTHON'

import gtfparse
import polars as pl

# Read GTF file
df_gtf = pl.from_pandas(gtfparse.read_gtf("${gtf}"))

# Filter for exon features only
exons = df_gtf.filter(pl.col("feature") == "exon")

# Convert to polars DataFrame for better performance
exons_pl = pl.DataFrame(exons)

# Add 0-based start for BED format
exons_pl = exons_pl.with_columns(
    start_0based=pl.col("start") - 1
)

# Number exons within each gene
exons_pl = (
    exons_pl
    .sort(["gene_id", "start"])
    .with_columns([
        pl.col("gene_id").cum_count().over("gene_id").alias("exon_number")
    ])
)

# Create name column in the required format
exons_pl = exons_pl.with_columns([
    (
        pl.col("gene_id") + "_exon_" +
            pl.col("exon_number").cast(pl.String) + "_0_chr" +
        pl.col("seqname") + "_" +
            (pl.col("start")).cast(pl.String) + "_" +
        pl.when(pl.col("strand") == "+").then(pl.lit("f")).otherwise(pl.lit("r"))
    ).alias("name")
])

# Select and rename columns for BED format
bed_df = exons_pl.select([
    pl.col("seqname").alias("chrom"),
    pl.col("start_0based").alias("start"),
    pl.col("end"),
    pl.col("name"),
    pl.lit(0).alias("score"),
    pl.col("strand")
])

# Write to BED file
bed_df.write_csv("${prefix}.bed", separator="\\t", include_header=False)
PYTHON
