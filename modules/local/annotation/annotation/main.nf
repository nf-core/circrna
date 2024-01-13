process ANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(intersection)
    tuple val(meta2), path(inv_intersection)
    
    output:
    tuple val(meta), path("annotation.bed"), emit: bed
    tuple val(meta2), path("annotation.tsv"), emit: locstring
    
    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import os

    columns = {
        0: 'chr',
        1: 'start',
        2: 'end',
        3: 'strand',
        12: 'attributes',
        13: 'overlap'
    }

    attributes = ['gene_id']

    df = pd.read_csv('$intersection', sep='\\t', header=None, usecols=columns.keys())
    df = df.rename(columns=columns)

    # Convert attributes to a dictionary
    df['attributes'] = df['attributes'].apply(lambda row: dict([[value.strip(r'\"') for value in entry.strip().split(' ')] for entry in row.split(';') if entry]))
    # Keep only the attributes we want
    df['attributes'] = df['attributes'].apply(lambda row: {key: row[key] for key in attributes})
    # Convert attributes to columns
    df = pd.concat([df.drop(['attributes'], axis=1), df['attributes'].apply(pd.Series)], axis=1)

    df = df.groupby(['chr', 'start', 'end', 'strand']).aggregate(lambda x: set(x))

    df['length'] = df.index.map(lambda row: row[2] - row[1])
    # Check if there is an overlap entry with the same length as the feature
    df['perfect'] = df.apply(lambda row: row['length'] in row['overlap'], axis=1)
    # Drop the overlap and length columns
    df = df.drop(['overlap', 'length'], axis=1)
    # Kee only unique gene_ids
    df['gene_id'] = df['gene_id'].apply(lambda row: ','.join(row))

    # If 'perfect': circRNA, else EIciRNA
    df['type'] = df['perfect'].apply(lambda perfect: 'circRNA' if perfect else 'EIciRNA')
    # Drop the perfect column
    df = df.drop(['perfect'], axis=1)

    if not os.stat('$inv_intersection').st_size == 0:
        df_nooverlap = pd.read_csv('$inv_intersection', sep='\\t', header=None)
        df_nooverlap.columns = ['chr', 'start', 'end', 'strand']
        df_nooverlap['type'] = 'ciRNA'
        df_nooverlap['gene_id'] = 'NA'

        # Set chr, start, end and strand as index
        df_nooverlap = df_nooverlap.set_index(['chr', 'start', 'end', 'strand'])

        df = pd.concat([df, df_nooverlap])

    df = df.sort_index()

    df.to_csv('annotation.bed', sep='\\t', index=True, header=False)

    df.index = df.index.map(lambda row: "{}:{}-{}:{}".format(*row))
    df.to_csv('annotation.tsv', sep='\\t', index=True, header=True)
    """
}