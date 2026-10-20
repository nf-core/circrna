process CIRCTOOLS_REFERENCE {
    tag "${genome}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f523ec9d03407da839bae40c76b9e49725d8b5273e9e167f7b5e081392f690d/data'
        : 'community.wave.seqera.io/library/circtools:2.0.5.post2--20ffd704d69916e2'}"

    input:
    val(genome)

    output:
    path "reference", emit: reference
    tuple val("${task.process}"), val('circtools'), eval('circtools -V'), emit: versions_circtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p temp
    export TMPDIR=./temp
    circtools nanopore -d -R reference/ -C ${genome}

    """

    stub:
    """
    mkdir -p reference

    """
}
