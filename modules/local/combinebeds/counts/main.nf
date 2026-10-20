process COMBINEBEDS_COUNTS {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b584f817e89b571a59ff2e6c887863e005c4ba32057e188ef2868942dc9a9a4/data'
        : 'community.wave.seqera.io/library/python_polars:a7975bec8c9bd1e2'}"

    input:
    tuple val(meta), val(aggregation), path(candidates), path(beds)
    val max_shift
    val consider_strand

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: combined, optional: true

    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template("counts.py")

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python --version 2>&1 | cut -d " " -f 2)"
        polars: "\$(python -c \'import polars; print(polars.__version__)\')"
    END_VERSIONS
    """
}
