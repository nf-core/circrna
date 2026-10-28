process CIRCTOOLS_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/51/51f47a2e827b3bc883ac39107f50e25e5ad8f014138300106751c6a444dbff3d/data'
        : 'community.wave.seqera.io/library/gtfparse_polars:58e453d2bfabb473'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('annotation.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: "\$(python3 -c 'import platform; print(platform.python_version())')"
        gtfparse: "\$(python3 -c 'import gtfparse; print(gtfparse.__version__)')"
        polars: "\$(python3 -c 'import polars; print(polars.__version__)')"
    END_VERSIONS
    """
}
