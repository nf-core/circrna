process CIRI_CIRIAS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f945f0dd8c450d3c034f23db6680a7c44e1447f790ed1cf58727cc9f1ed02e72/data'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--eca1ab13065c3c7d'}"

    input:
    tuple val(meta), path(ciri2), path(sam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}_jav.list"), emit: list
    tuple val(meta), path("${prefix}_library_length.list"), emit: library_length
    tuple val("${task.process}"), val('ciri-as'), val('1.2'), emit: versions_ciri_as, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI-AS -S ${sam} -C ${ciri2} -F ${fasta} -A ${gtf} -O ${prefix} -D yes ${args}

    # The library length file stays even if it is empty: CIRI-vis needs it via -l to write the circRNA sequences

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_jav.list
    touch ${prefix}_library_length.list

    """
}
