process UNIFY {
    tag "${meta.id}"
    label 'process_single'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/polars:1.24.0--800cd3e4ff805434'
        : 'community.wave.seqera.io/library/polars:1.24.0--2d2d323e8514e707'}"

    input:
    tuple val(meta), path(reads), path(coordinates), path(counts)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: bed
    path "versions.yml"        , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "bed"
    template 'unify.py'
}
