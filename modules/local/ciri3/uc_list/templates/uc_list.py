#!/usr/bin/env python3

import platform

import polars as pl
import yaml


def main(bed_path, output_path):
    schema = {
        "chrom": pl.String, "start": pl.Int64, "end": pl.Int64,
        "name": pl.String, "score": pl.String, "strand": pl.String,
        "thick_start": pl.String, "thick_end": pl.String, "item_rgb": pl.String,
        "block_count": pl.String, "block_sizes": pl.String, "block_starts": pl.String,
    }
    try:
        bed = pl.read_csv(
            bed_path, separator="\\t", has_header=False, schema=schema,
            comment_prefix="#", quote_char=None, truncate_ragged_lines=True,
            raise_if_empty=False,
        )
    except pl.exceptions.ComputeError as error:
        raise ValueError("BED12 coordinates must be integers") from error
    if bed.is_empty():
        raise ValueError("consensus BED12 contains no circRNA candidates")
    if bed["block_starts"].is_null().any():
        raise ValueError("BED12 input: expected at least 12 columns")
    candidates = bed.select("chrom", (pl.col("start") + 1).alias("start"), "end")
    if candidates.null_count().sum_horizontal().item() or candidates.filter(
        (pl.col("chrom") == "") | (pl.col("start") < 1) | (pl.col("end") < pl.col("start"))
    ).height:
        raise ValueError("BED12 input: invalid interval")
    if candidates.is_duplicated().any():
        raise ValueError("BED12 input: duplicate consensus interval")
    candidates.write_csv(output_path, separator="\\t", include_header=False)


if __name__ == "__main__":
    main("${bed}", "${prefix}.ciri3.list")
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "polars": pl.__version__,
            "pyyaml": yaml.__version__,
        }
    }
    with open("versions.yml", "w", encoding="utf-8") as version_file:
        yaml.safe_dump(versions, version_file)
