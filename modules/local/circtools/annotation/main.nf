process CIRCTOOLS_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/gtfparse_pyyaml:e3b66ef529a096e7'
        : 'community.wave.seqera.io/library/gtfparse_pyyaml:508c7607beb17f69'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml"        , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'annotation.py'
}
