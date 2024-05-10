process MERGE_EXPERIMENTS {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.62.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-rtracklayer:1.62.0--r43ha9d7317_0' }"

    input:
    tuple val(meta),  path(experiments)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(tpm)

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
