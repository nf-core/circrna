process FISHPOND_SWISH {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-fishpond:2.8.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-fishpond:2.8.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(experiment)
    val(column)

    output:
    tuple val(meta), path("${meta.id}.swish.rds"), emit: swish
    path "versions.yml"                    , emit: versions

    script:
    template "swish.r"

    stub:
    """
    touch ${meta.id}.RDS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-fishpond: \$(Rscript -e "library(fishpond); cat(as.character(packageVersion('fishpond')))")
    END_VERSIONS
    """
}
