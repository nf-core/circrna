process MAPSPLICE {
    tag "$meta.id"
    label 'process_high'

/*         conda (params.enable_conda ? "bioconda::mapsplice=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }" */

}
