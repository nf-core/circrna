#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PYTHON'

from io import StringIO

import polars as pl


def main(bed_path, output_path):
    schema = {
        "chrom": pl.String, "start": pl.Int64, "end": pl.Int64,
        "name": pl.String, "score": pl.String, "strand": pl.String,
        "thick_start": pl.String, "thick_end": pl.String, "item_rgb": pl.String,
        "block_count": pl.String, "block_sizes": pl.String, "block_starts": pl.String,
    }
    # Polars 2 infers the CSV width from the first record, including blank lines.
    # Normalize BED comments/blank lines before typed parsing, retaining row order.
    with open(bed_path, encoding="utf-8") as handle:
        records = StringIO("".join(line for line in handle if line.strip() and not line.lstrip().startswith("#")))
    try:
        bed = pl.read_csv(
            records, separator="\\t", has_header=False, schema=schema,
            comment_prefix="#", quote_char=None, truncate_ragged_lines=True,
            raise_if_empty=False, extra_columns="ignore",
        )
    except pl.exceptions.SchemaError as error:
        raise ValueError("BED12 input: expected at least 12 columns") from error
    except pl.exceptions.ComputeError as error:
        raise ValueError("BED12 coordinates must be integers") from error
    bed = bed.filter(~pl.all_horizontal(pl.all().is_null()))
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
PYTHON
