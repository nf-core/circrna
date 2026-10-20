process CIRCMINER_ORDERGTF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/47c5d77ca8b261780a14c2f3ef37ff3662bdec43474e9675d791cfb8ce383d55/data'
        : 'community.wave.seqera.io/library/python_gffutils:65d4e6d2f4b119ea'}"

    input:
    tuple val(meta), path(gtf, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.circminer"
    template('ordergtf.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.circminer"
    """
    touch ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: "\$(python3 -c 'import platform; print(platform.python_version())')"
        gffutils: "\$(python3 -c 'import gffutils; print(gffutils.__version__)')"
    END_VERSIONS
    """
}
