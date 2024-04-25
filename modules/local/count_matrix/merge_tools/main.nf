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

    path  "versions.yml"                          , emit: versions

    script:
    """
    merge_tools.py --beds ${beds} --tool_filter ${tool_filter} --duplicates_fun ${duplicates_fun} --output ${meta.id}.merged.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
