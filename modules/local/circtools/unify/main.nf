process UNIFY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc8e5ab4ba32d58b69d5395ae4b39f0fde2467e5a81c8a6895d664309bb27498/data'
        : 'community.wave.seqera.io/library/polars:2.0.0--2895a858bbcef9a5'}"

    input:
    tuple val(meta), path(reads), path(coordinates), path(counts)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('unify.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: "\$(python3 -c 'import platform; print(platform.python_version())')"
        polars: "\$(python3 -c 'import polars; print(polars.__version__)')"
    END_VERSIONS
    """
}
