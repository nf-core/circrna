process CIRI2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f945f0dd8c450d3c034f23db6680a7c44e1447f790ed1cf58727cc9f1ed02e72/data'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--eca1ab13065c3c7d'}"

    input:
    tuple val(meta), path(sam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('ciri2'), val('2.0.6'), emit: versions_ciri2, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI -I ${sam} -O ${prefix}.txt -F ${fasta} -A ${gtf} -T ${task.cpus} -0 ${args}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    """
}
