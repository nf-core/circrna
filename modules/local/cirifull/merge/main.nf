process CIRIFULL_MERGE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f945f0dd8c450d3c034f23db6680a7c44e1447f790ed1cf58727cc9f1ed02e72/data'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--eca1ab13065c3c7d'}"

    input:
    tuple val(meta), path(ciri2), path(ciri_as), path(ro2)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}_merge_circRNA_detail.anno"), emit: anno
    tuple val("${task.process}"), val('ciri-full'), val('2.1.2'), emit: versions_ciri_full, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI-full Merge -a ${gtf} -r ${fasta} -c ${ciri2} -as ${ciri_as} -ro ${ro2} -o ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merge_circRNA_detail.anno
    """
}
