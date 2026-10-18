#!/usr/bin/env python3

import platform
import yaml

circ_list_path = "${circ_list}"
bsj_path = "${bsj_matrix}"
fsj_path = "${fsj_matrix}"
prefix = "${prefix}"


def fail(message):
    raise ValueError(message)


def parse_candidates(path):
    candidates = []
    seen = set()
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) != 3 or not fields[0]:
                fail(f"candidate list row {line_number}: expected chrom, start, end")
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"candidate list row {line_number}: coordinates must be integers") from error
            if start < 1 or end <= start:
                fail(f"candidate list row {line_number}: invalid 1-based interval {chrom}:{start}-{end}")
            key = (chrom, start, end)
            if key in seen:
                fail(f"duplicate candidate at row {line_number}: {key}")
            seen.add(key)
            candidates.append(key)
    if not candidates:
        fail("candidate list is empty")
    return candidates, seen


def parse_matrix(path, expected_header, candidates):
    counts = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\\r\\n").split("\\t")
        if header != ["circRNA", expected_header]:
            fail(f"{path}: expected header circRNA\\t{expected_header}, got {header!r}")
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\\r\\n").split("\\t")
            if len(fields) != 4:
                fail(f"{path} row {line_number}: expected chrom, start, end, count")
            chrom = fields[0]
            try:
                key = (chrom, int(fields[1]), int(fields[2]))
                count = int(fields[3])
            except ValueError as error:
                raise ValueError(f"{path} row {line_number}: coordinates/count must be integers") from error
            if key not in candidates:
                fail(f"{path} row {line_number}: CIRI3 returned a coordinate outside the requested consensus: {key}")
            if count < 0:
                fail(f"{path} row {line_number}: count must not be negative")
            if key in counts:
                fail(f"{path} row {line_number}: duplicate matrix coordinate {key}")
            counts[key] = count
    return counts


candidates, candidate_set = parse_candidates(circ_list_path)
bsj_counts = parse_matrix(bsj_path, "BSJ", candidate_set)
fsj_counts = parse_matrix(fsj_path, "FSJ", candidate_set)

for suffix, counts in (("bsj", bsj_counts), ("fsj", fsj_counts)):
    with open(f"{prefix}.{suffix}.tsv", "w", encoding="utf-8") as output:
        output.write("circRNA_ID\\tcount\\n")
        for chrom, start, end in candidates:
            key = (chrom, start, end)
            output.write(f"{chrom}:{start}|{end}\\t{counts.get(key, 0)}\\n")

versions = {
    "${task.process}": {
        "python": platform.python_version(),
    }
}
with open("versions.yml", "w", encoding="utf-8") as version_file:
    yaml.safe_dump(versions, version_file)
