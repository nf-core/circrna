process EXTRACT_TRANSCRIPTOME {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.transcriptome.bed"), emit: bed

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

    df = pd.read_csv("${gtf}", sep="\\t", comment="#", names=columns, header=None)

    # Keep only transcript features
    df = df[df['feature'] == 'transcript']

    # Convert attributes to a dictionary
    df['attributes'] = df['attributes'].apply(lambda row: dict([[value.strip(r'"') for value in entry.strip().split(' ')] for entry in row.split(';') if entry]))
    
    # Extract gene_id and transcript_id
    df['gene_id'] = df['attributes'].apply(lambda row: row['gene_id'])
    df['transcript_id'] = df['attributes'].apply(lambda row: row['transcript_id'])
    df['name'] = df['gene_id'] + ":" + df['transcript_id']

    # Format as bed
    df = df[['seqname', 'start', 'end', 'name', 'score', 'strand']]

    df.to_csv("${meta.id}.transcriptome.bed", sep="\\t", index=False, header=False)
    """
}
