process CIRCTOOLS_NANOPORE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b2d16e12d96fbed5965f922347a73c59232ccf4420325577e6c1a20e788f0a31/data'
        : 'community.wave.seqera.io/library/bedtools_circtools_nanofilt_pblat_samtools:7a25c5c7268bf53b'}"

    input:
    tuple val(meta), path(reads)
    path(reference)
    val(genome)

    output:
    tuple val("${task.process}"), val('circtools'), eval('circtools -V'), emit: versions_circtools, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p temp
    export TMPDIR="\$(pwd)/temp"

    circtools nanopore \
        -r \
        -s ${reads} \
        -R ${reference} \
        -C ${genome} \
        -O ${prefix} \
        -t ${task.cpus}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    """
}
