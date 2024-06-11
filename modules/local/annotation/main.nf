process ANNOTATION {
    tag "$meta.id:$meta.tool"
    label 'process_single'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(intersection)
    val(exon_boundary)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}.gtf"), emit: gtf

    path "versions.yml"                   , emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'annotation.py'
}
