process PNG_JSON {
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(png)
    val(title)
    val(description)

    output:
    path("*_mqc.json")  , emit: report
    path("versions.yml"), emit: versions

    script:
    template "json.py"
}
