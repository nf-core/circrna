process BENCHMARKING_MULTIQC {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(jaccard)

    output:
    path("*_mqc.json")  , emit: report
    path("versions.yml"), emit: versions

    script:
    template "benchmarking.py"
}
