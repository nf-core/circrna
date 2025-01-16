#!/usr/bin/env python3
import platform
import polars as pl
import yaml

# paths = "${bindingsites}".split(" ")
paths = []
bindingsites = "${bindingsites}"

j = 0
for i in range(len(bindingsites)):
    if bindingsites[i] == " " or i == len(bindingsites) - 1:
        paths.append(bindingsites[j : i + 1])
        j = i + 1

print(paths[1:3])
print(len(paths))

df = pl.scan_csv(paths,
                 separator="\\t",
                 has_header=False,
                 new_columns=['mirna', 'target', 'start', 'end', 'tool'])

df = df.select(["mirna", "target", "tool"])

df = df.group_by(['mirna', 'target']).agg(pl.col("tool").n_unique())

df = df.filter(pl.col("tool") >= int("${min_tools}")) \
    .select(["mirna", "target"])

df = df.collect()

df.write_csv('${meta.id}.majority.tsv', separator='\\t', include_header=False)

# Create targets file

df = df.group_by('mirna').agg(pl.col("target").str.concat(","))

df.write_csv('${meta.id}.targets.tsv', separator='\\t', include_header=False)

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "polars": pl.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))

