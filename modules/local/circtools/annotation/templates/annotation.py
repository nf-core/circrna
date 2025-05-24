#! /usr/bin/env python3

import platform

import gtfparse
import polars as pl
import yaml

# Read GTF file
df_gtf = gtfparse.read_gtf("${gtf}")

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
        pl.col("gene_id").cumcount().over("gene_id").alias("exon_number")
    ])
)

# Create name column in the required format
exons_pl = exons_pl.with_columns([
    (
        pl.col("gene_id") + "_exon_" +
        pl.col("exon_number").cast(pl.Utf8) + "_0_chr" +
        pl.col("seqname") + "_" +
        (pl.col("start")).cast(pl.Utf8) + "_" +
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

# Log versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
