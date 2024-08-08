process COMPUTE_CORRELATIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-fishpond:2.8.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-fishpond:2.8.0--r43hdfd78af_0' }"

    input:
    tuple val(meta),  path(bindingsites)
    tuple val(meta2), path(mirna_expression)
    tuple val(meta3), path(transcript_rds)

    output:
    tuple val(meta), path("*.tsv"), emit: correlations, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'compute_correlations.R'

    stub:
    """
    touch ${meta.id}.circrna_correlation.tsv
    """
}
