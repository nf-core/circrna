#! /usr/bin/env python3

import platform

import polars as pl

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

df_counts = pl.scan_csv("${counts}", separator="\\t")
df_counts = df_counts.rename({"paired.junctions": "Count"})

join_columns = ['Chr', 'Start', 'End']
df_coordinates = pl.scan_csv("${coordinates}", separator="\\t")
df_coordinates = df_coordinates.select(join_columns + ['Strand'])

df_counts = df_counts.join(df_coordinates, on=join_columns, how='left')

df_reads = pl.scan_csv("${reads}", separator="\\t", has_header=False, schema = {
    "chr_donorA": pl.Utf8,              # Chromosome of the donor site
    "start1": pl.Int64,                 # Genomic position of the donor site
    "strand1": pl.Utf8,                 # Strand of the donor site (+ or -)
    "chr_acceptorB": pl.Utf8,           # Chromosome of the acceptor site
    "start2": pl.Int64,                 # Genomic position of the acceptor site
    "strand2": pl.Utf8,                 # Strand of the acceptor site (+ or -)
    "junction_type": pl.Int64,          # Type of junction (0-6)
    "repeat_left_junction": pl.Int64,   # Bases extending left of junction
    "repeat_right_junction": pl.Int64,  # Bases extending right of junction
    "read_id": pl.Utf8,                 # Name of the read
    "first_base_read1": pl.Int64,       # Position of first base of read1
    "CIGAR_read1": pl.Utf8,             # CIGAR string of read1
    "first_base_read2": pl.Int64,       # Position of first base of read2
    "CIGAR_read2": pl.Utf8,             # CIGAR string of read2
    "sample": pl.Utf8
})
df_reads = df_reads.filter((pl.col("chr_donorA") == pl.col("chr_acceptorB")) & (pl.col("strand1") == pl.col("strand2")))
df_reads = df_reads.with_columns(
    Chr=pl.col("chr_donorA"),
    Strand=pl.col("strand1"),
    # Start is first_base_read1 if strand is +, otherwise it is first_base_read2
    Start=pl.when(pl.col("strand1") == "+").then(pl.col("first_base_read2")).otherwise(pl.col("first_base_read1")),
    # End is first_base_read2 if strand is +, otherwise it is first_base_read1
    End=pl.when(pl.col("strand1") == "+").then(pl.col("start1")).otherwise(pl.col("start2"))-1,
)
df_reads = df_reads.drop(["chr_donorA", "start1", "strand1", "chr_acceptorB", "start2", "strand2"])
df_counts = df_counts.join(df_reads, on=["Chr", "Start", "End"], how="left")
df_counts = df_counts.group_by(["Chr", "Start", "End", "Count", "Strand"]).agg(
    Reads=pl.col("read_id").unique().str.join(",")
)
df_counts = df_counts.with_columns(
    Start=pl.col("Start")-1
)

df_counts = df_counts.collect()
df_counts = df_counts.with_columns(
    Name=pl.lit("FUSIONJUNC_") + pl.int_range(0, len(df_counts)).cast(pl.Utf8) + pl.lit("/") + pl.col("Count").cast(pl.Utf8)
)
df_counts = df_counts.select(["Chr", "Start", "End", "Name", "Count", "Strand", "Reads"])
df_counts.write_csv("${prefix}.${suffix}", separator="\\t", include_header=False)

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
