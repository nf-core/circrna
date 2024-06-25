process AVERGAGE_TSV {
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"


    input:
    path(tsv)

    output:
    path(tsv)

    script:
    template "average.py"
}
