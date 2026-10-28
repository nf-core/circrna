#!/usr/bin/env python3

import platform
import polars as pl

length_files = "${lengths}".split()
df = pl.scan_csv(length_files, separator='\\t', has_header=False, new_columns=["read_id", "length"])

# Find the 5th percentile
percentile = (df.select('length')
                .quantile(0.05, interpolation='lower')
                .collect()['length'][0])

# Write the result to a file
with open('${prefix}.txt', 'w') as f:
    f.write(str(int(percentile)))
with open("versions.yml", "w", encoding="utf-8") as handle:
    print("${task.process}:", file=handle)
    print(f'    python: "{platform.python_version()}"', file=handle)
    print(f'    polars: "{pl.__version__}"', file=handle)
