process ANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(intersection)
    val(exon_boundary)

    output:
    tuple val(meta), path("${meta.id}.annotation.bed"), emit: bed

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import os

    columns = {
        0: 'chr',
        1: 'start',
        2: 'end',
        3: 'name',
        4: 'score',
        5: 'strand',
        9: 'tx_start',
        10: 'tx_end',
        14: 'attributes'
    }

    attributes = ['gene_id', 'transcript_id']

    exon_boundary = ${exon_boundary}

    df = pd.read_csv("${intersection}", sep="\\t", header=None, usecols=columns.keys())
    df = df.rename(columns=columns)

    # Extract circRNAs without match
    mask = df['tx_start'] == -1
    df_nomatch = df[mask]
    df = df[~mask]
    df_nomatch['type'] = 'intergenic-circRNA'
    df_nomatch['gene_id'] = 'NaN'
    df_nomatch['transcript_id'] = 'NaN'

    # Convert attributes to a dictionary
    df['attributes'] = df['attributes'].apply(lambda row: dict([[value.strip(r'"') for value in entry.strip().split(' ')] for entry in row.split(';') if entry]))
    # Keep only the attributes we want
    df['attributes'] = df['attributes'].apply(lambda row: {key: row[key] for key in attributes})
    # Convert attributes to columns
    df = pd.concat([df.drop(['attributes'], axis=1), df['attributes'].apply(pd.Series)], axis=1)

    df['any_outside'] = (df['start'] < df['tx_start'] - exon_boundary) | (df['end'] > df['tx_end'] + exon_boundary)
    # Perfect is inverse of any_outside
    df['perfect'] = ~df['any_outside']
    # Drop any_outside
    df = df.drop(['any_outside', 'tx_start', 'tx_end'], axis=1)

    df = df.groupby(['chr', 'start', 'end', 'strand']).aggregate({
        'name': lambda x: x.iloc[0],
        'score': lambda x: x.iloc[0],
        'gene_id': lambda x: list(x),
        'transcript_id': lambda x: list(x),
        'perfect': lambda x: list(x)
    })

    def filter_perfect(row, col):
        if any(row['perfect']):
            matching_values = [value for value, perfectness in zip(row[col], row['perfect']) if perfectness]
        else:
            matching_values = row[col]
        return ",".join(set(matching_values))

    df['type'] = df['perfect'].apply(lambda x: 'circRNA' if any(x) else 'EI-circRNA')
    df['gene_id'] = df.apply(lambda row: filter_perfect(row, 'gene_id'), axis=1)
    df['transcript_id'] = df.apply(lambda row: filter_perfect(row, 'transcript_id'), axis=1)
    # Drop perfect
    df = df.drop(['perfect'], axis=1)

    df = df.reset_index()
    df_nomatch = df_nomatch.reset_index()
    bed_order = ['chr', 'start', 'end', 'name', 'score', 'strand', 'type', 'gene_id', 'transcript_id']
    df = df[bed_order]
    df_nomatch = df_nomatch[bed_order]

    df = pd.concat([df, df_nomatch], axis=0)

    # Sort by chr, start, end
    df = df.sort_values(['chr', 'start', 'end'])

    df.to_csv('${meta.id}.annotation.bed', sep='\\t', index=False, header=False)
    """
}
