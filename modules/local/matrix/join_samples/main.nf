process JOIN_SAMPLES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0a88e97bf9947f619334029862edf5b0a560196105f003aa9d0d5905f58d2a40/data'
        : 'community.wave.seqera.io/library/join_samples:f5ee4ff1280edbba'}"

    input:
    tuple val(meta), val(samples), path(matrices)
    val metacols

    output:
    tuple val(meta), path("${meta.id}.joined.tsv"), emit: joined

    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('join.py')

    stub:
    """
    touch ${meta.id}.joined.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python --version 2>&1 | cut -d " " -f 2)"
        polars: "\$(python -c \'import polars; print(polars.__version__)\')"
    END_VERSIONS
    """
}
