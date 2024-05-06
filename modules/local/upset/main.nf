process UPSET {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3 conda-forge::numpy=1.20.* conda-forge::pandas=1.2.* conda-forge::upsetplot=0.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"
    input:
    tuple val(meta), val(tools), path(beds)

    when:
    task.ext.when == null || task.ext.when

    output:
    tuple val(meta), path("*.png"), emit: plot
    path "*.upset_mqc.json"       , emit: multiqc
    path "versions.yml"           , emit: versions

    script:
    template "upset.py"
}
