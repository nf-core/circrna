process MERGE_EXPERIMENTS {
    tag "${meta.id}"
    label "process_medium"
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/ee16c4b13326f90cb63ced65618276af1f4dd061f0d405ce87d9c78c863e4c75/data'
        : 'community.wave.seqera.io/library/merge_experiments:85f748798d2b090d'}"

    input:
    tuple val(meta), path(experiments)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(tpm)

    output:
    tuple val(meta), path("${meta.id}.merged.rds"), emit: merged
    tuple val("${task.process}"), val('bioconductor-summarizedexperiment'), eval('Rscript -e "cat(as.character(packageVersion(\'SummarizedExperiment\')))"'), emit: versions_bioconductor_summarizedexperiment, topic: versions
    tuple val("${task.process}"), val('bioconductor-rtracklayer'), eval('Rscript -e "cat(as.character(packageVersion(\'rtracklayer\')))"'), emit: versions_bioconductor_rtracklayer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("merge_experiments.sh")

    stub:
    """
    touch ${meta.id}.merged.rds

    """
}
