process MERGE_TOOLS {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(beds)
    val(tool_filter)
    val(duplicates_fun)

    output:
    tuple val(meta), path("${meta.id}.merged.bed"), emit: merged

    path "versions.yml"                           , emit: versions

    script:
    template 'merge_tools.py'
}
