process CIRIFULL_RO2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f945f0dd8c450d3c034f23db6680a7c44e1447f790ed1cf58727cc9f1ed02e72/data'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--eca1ab13065c3c7d'}"

    input:
    tuple val(meta), path(sam), val(length)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}_ro2_info.list"), emit: list
    tuple val("${task.process}"), val('ciri-full'), val('2.1.2'), emit: versions_ciri_full, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI-full RO2 -r ${fasta} -s ${sam} -l ${length} -o ${prefix} $args

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ro2_info.list

    """
}
