process MERGE_EXPERIMENTS {
    tag "${meta.id}"
    label "process_medium"

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
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("merge_experiments.r")

    stub:
    """
    touch ${meta.id}.merged.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: "\$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")"
        bioconductor-rtracklayer: "\$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")"
    END_VERSIONS

    """
}
