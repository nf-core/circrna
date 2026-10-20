process PSIRC_BSJ_UNIFY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f26bc70bd7dfff0a45e25dec0bd8f1a3dfb146aab84b68437aac652f74b7888/data'
        : 'community.wave.seqera.io/library/psirc_bsj_unify:a8c5b7b17067ac40'}"

    input:
    tuple val(meta), path(bed, stageAs: 'input/*')
    tuple val(transcriptome_meta), path(transcriptome)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval('python --version | cut -d " " -f2'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template('psirc_bsj_unify.sh')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.bed
    """
}
