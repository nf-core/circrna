process FIND_CIRC_ANCHORS {
    tag "${meta.id}"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f572ecb116eee7e428c2f34e11a9c98d00b922893b1e97617a4875743ab967e/data'
        : 'community.wave.seqera.io/library/find_circ_anchors:e16458c7823af09d'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}_anchors.qfa.gz"), emit: anchors
    tuple val("${task.process}"), val('find_circ'), val('1.2'), emit: versions_find_circ, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    unmapped2anchors.py ${bam} | gzip > ${prefix}_anchors.qfa.gz

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_anchors.qfa.gz

    """
}
