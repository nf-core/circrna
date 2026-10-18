#!/usr/bin/env python3

import platform

import polars as pl
import yaml

circ_list_path = "${circ_list}"
bsj_path = "${bsj_matrix}"
fsj_path = "${fsj_matrix}"
prefix = "${prefix}"


COORDINATES = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64}


def parse_candidates(path):
    try:
        candidates = pl.read_csv(
            path, separator="\\t", has_header=False, schema=COORDINATES,
            quote_char=None, raise_if_empty=False,
        )
    except pl.exceptions.ComputeError as error:
        raise ValueError("candidate list: expected chrom, start, end with integer coordinates") from error
    if candidates.is_empty():
        raise ValueError("candidate list is empty")
    if candidates.null_count().sum_horizontal().item() or candidates.filter(
        (pl.col("chrom") == "") | (pl.col("start") < 1) | (pl.col("end") <= pl.col("start"))
    ).height:
        raise ValueError("candidate list: invalid 1-based interval")
    if candidates.is_duplicated().any():
        raise ValueError("duplicate candidate in consensus list")
    return candidates


def parse_matrix(path, expected_header, candidates):
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\\r\\n").split("\\t")
    if header != ["circRNA", expected_header]:
        raise ValueError(f"{path}: expected header circRNA\\t{expected_header}, got {header!r}")
    try:
        counts = pl.read_csv(
            path, separator="\\t", has_header=False, skip_rows=1,
            schema={**COORDINATES, "count": pl.Int64}, quote_char=None,
            raise_if_empty=False,
        )
    except pl.exceptions.ComputeError as error:
        raise ValueError(f"{path}: coordinates/count must be integers") from error
    if counts.null_count().sum_horizontal().item():
        raise ValueError(f"{path}: expected chrom, start, end, count without missing values")
    if counts.filter(pl.col("count") < 0).height:
        raise ValueError(f"{path}: count must not be negative")
    if counts.select(*COORDINATES).is_duplicated().any():
        raise ValueError(f"{path}: duplicate matrix coordinate")
    if counts.join(candidates, on=list(COORDINATES), how="anti").height:
        raise ValueError(f"{path}: CIRI3 returned a coordinate outside the requested consensus")
    return counts


candidates = parse_candidates(circ_list_path)

for suffix, path in (("bsj", bsj_path), ("fsj", fsj_path)):
    counts = parse_matrix(path, suffix.upper(), candidates)
    candidates.join(counts, on=list(COORDINATES), how="left", maintain_order="left").select(
        pl.format("{}:{}|{}", "chrom", "start", "end").alias("circRNA_ID"),
        pl.col("count").fill_null(0),
    ).write_csv(f"{prefix}.{suffix}.tsv", separator="\\t")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
        "pyyaml": yaml.__version__,
    }
}
with open("versions.yml", "w", encoding="utf-8") as version_file:
    yaml.safe_dump(versions, version_file)
