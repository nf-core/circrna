#!/usr/bin/env python

import platform
import yaml

import polars as pl

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "polars": pl.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))

# Main

exons_only = bool("${exons_only}")

columns = ['chr', 'start', 'end', 'name', 'score', 'strand',
            'thickStart', 'thickEnd', 'itemRgb',
            'exonCount', 'exonSizes', 'exonStarts',
            'readNumber', 'circType', 'gene', 'transcript',
            'index', 'flankIntron'
            ]
try:
    df = pl.scan_csv('${bed12}', separator='\\t', has_header=False, new_columns=columns, raise_if_empty=True)
except pl.exceptions.NoDataError:
    with open('${prefix}.${suffix}', 'w') as f:
        f.write('')
    exit(0)

df = df.with_columns(
    attributes = pl.lit('gene_id "') + pl.col('gene') + pl.lit('"; transcript_id "') + pl.col('name') + pl.lit('";'),
    source = pl.lit('nf-core/circrna')
)

if exons_only:
    df_exons = df.with_columns(
        exonSizes = pl.col('exonSizes').str.split(','),
        exonStarts = pl.col('exonStarts').str.split(',')
    ).explode('exonSizes', 'exonStarts')
    df_exons = df_exons.with_columns(
        type = pl.lit('exon'),
        phase = pl.lit('.'),
        start = pl.col('start') + pl.col('exonStarts').cast(int)
    ).with_columns(
        end = pl.col('start') + pl.col('exonSizes').cast(int)
    )
    df_exons = df_exons.select(
        'chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'
    )
    df_cds = df_exons.clone()
    df_cds = df_cds.with_columns(
        type = pl.lit('CDS'),
        phase = pl.lit('.')
    )
else:
    df_exons = df.clone()
    df_exons = df_exons.with_columns(
        type = pl.lit('exon'),
        phase = pl.lit('.')
    )
    df_exons = df_exons.select(
        'chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'
    )
    df_cds = df_exons.clone()
    df_cds = df_cds.with_columns(
        type = pl.lit('CDS'),
        phase = pl.lit('.')
    )

df_combined = pl.concat([df_exons, df_cds])
df_combined = df_combined.sort('chr', 'start', 'end')

df_combined.collect().write_csv('${prefix}.${suffix}', separator='\\t', include_header=False, quote_style="never")
