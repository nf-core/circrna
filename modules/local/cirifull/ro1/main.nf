process CIRIFULL_RO1 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f945f0dd8c450d3c034f23db6680a7c44e1447f790ed1cf58727cc9f1ed02e72/data'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--eca1ab13065c3c7d'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_ro1_align.txt"), emit: align
    tuple val(meta), path("${prefix}_ro1.fq.gz"), emit: fastq
    tuple val("${task.process}"), val('ciri-full'), val('2.1.2'), emit: versions_ciri_full, topic: versions
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI-full RO1 -1 ${reads[0]} -2 ${reads[1]} -o ${prefix} $args
    gzip ${prefix}_ro1.fq

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ro1_align.txt
    echo "" | gzip > ${prefix}_ro1.fq.gz

    """
}
