#! /usr/bin/env python3

import platform
import polars as pl
import yaml

length_files = "${lengths}".split()
df = pl.scan_csv(length_files, separator='\\t', has_header=False, new_columns=["read_id", "length"])

# Find the 5th percentile
percentile = (df.select('length')
                .quantile(0.05, interpolation='lower')
                .collect()['length'][0])

# Write the result to a file
with open('${prefix}.txt', 'w') as f:
    f.write(str(int(percentile)))

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
