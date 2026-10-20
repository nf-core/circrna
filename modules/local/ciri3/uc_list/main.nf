process CIRI3_UC_LIST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc24a035d2b10ae5cddc290f972da775bed663c0a9f3d046386412b2922112ad/data'
        : 'community.wave.seqera.io/library/python_polars:708bbfa3590a8658'}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${prefix}.ciri3.list"), emit: list
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('uc_list.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ciri3.list
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: "\$(python3 -c 'import platform; print(platform.python_version())')"
        polars: "\$(python3 -c 'import polars; print(polars.__version__)')"
    END_VERSIONS
    """
}
