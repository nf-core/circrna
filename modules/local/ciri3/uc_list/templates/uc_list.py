#!/usr/bin/env python3

import platform
import yaml


def main(bed_path, output_path):
    seen = set()
    with open(bed_path, encoding="utf-8") as bed, open(output_path, "w", encoding="utf-8") as output:
        for line_number, line in enumerate(bed, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\\r\\n").split("\\t")
            if len(fields) < 12:
                raise ValueError(f"BED12 row {line_number}: expected at least 12 columns, found {len(fields)}")
            chrom = fields[0]
            try:
                start0, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"BED12 row {line_number}: coordinates must be integers") from error
            if not chrom or start0 < 0 or end <= start0:
                raise ValueError(f"BED12 row {line_number}: invalid interval {chrom}:{start0}-{end}")
            candidate = (chrom, start0 + 1, end)
            if candidate in seen:
                raise ValueError(f"BED12 row {line_number}: duplicate consensus interval {candidate}")
            seen.add(candidate)
            output.write("\\t".join(map(str, candidate)) + "\\n")
    if not seen:
        raise ValueError("consensus BED12 contains no circRNA candidates")


if __name__ == "__main__":
    main("${bed}", "${prefix}.ciri3.list")
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
        }
    }
    with open("versions.yml", "w", encoding="utf-8") as version_file:
        yaml.safe_dump(versions, version_file)
