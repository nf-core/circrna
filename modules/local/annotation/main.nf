process ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(gtf_intersection), path(db_intersections)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}.gtf"), emit: gtf

    path "versions.yml"                   , emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'annotation.py'
}
