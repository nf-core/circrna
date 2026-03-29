process TXIMETA_TXIMETA {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta%3A1.20.1--r43hdfd78af_1' :
        'biocontainers/bioconductor-tximeta:1.20.1--r43hdfd78af_1' }"

    input:
    tuple val(meta), path("quants/*")
    val quant_type

    output:
    tuple val(meta), path("*.rds"), emit: se
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'tximeta.r'

    stub:
    """
    touch ${meta.id}.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
