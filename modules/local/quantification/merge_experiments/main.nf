process MERGE_EXPERIMENTS {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(experiments)
    tuple val(meta2), path(phenotype)

    output:
    tuple val(meta), path("${meta.id}.merged.rds"), emit: merged
    path "versions.yml"                           , emit: versions

    script:
    template "merge_experiments.r"

    stub:
    """
    touch ${meta.id}.merged.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
