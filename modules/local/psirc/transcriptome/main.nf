process PSIRC_TRANSCRIPTOME {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5afa9a01488ef93da1c7d72dfb0f0c9d885e5933aa67aa4ccb52f3c98eeb0c/data'
        : 'community.wave.seqera.io/library/psirc:1.0.0--4c24a8f727283097'}"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: transcriptome
    tuple val("${task.process}"), val('psirc'), val('1.0.0'), emit: versions_psirc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    suffix = 'fasta'
    template('transcriptome.sh')

    stub:
    prefix = task.ext.prefix ?: meta.id
    suffix = 'fasta'
    """
    touch ${prefix}.${suffix}
    """
}
