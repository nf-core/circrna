process FIND_CIRC_ANCHORS {
    tag "$meta.id"
    label "process_high"

    container 'barryd237/find_circ'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}_anchors.qfa.gz"), emit: anchors
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    unmapped2anchors.py $bam | gzip > ${prefix}_anchors.qfa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_circ: $VERSION
    END_VERSIONS
    """
}
