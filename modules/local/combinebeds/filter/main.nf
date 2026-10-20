process COMBINEBEDS_FILTER {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/48af262bc66e75295b1a314cbf45b01bded343911ffeaa0357009f7370cb4573/data'
        : 'community.wave.seqera.io/library/python_polars_upsetplot_pandas_matplotlib:a16ba3398a0ae7a3'}"

    input:
    tuple val(meta), path(beds)
    val max_shift
    val consider_strand
    val min_tools
    val min_samples
    val suffix

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: combined, optional: true
    path "*.png", emit: plots, optional: true
    path "*.json", emit: multiqc, optional: true

    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template("filter.py")

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${suffix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python --version 2>&1 | cut -d " " -f 2)"
        polars: "\$(python -c \'import polars; print(polars.__version__)\')"
        upsetplot: "\$(python -c \'import upsetplot; print(upsetplot.__version__)\')"
        matplotlib: "\$(python -c \'import matplotlib; print(matplotlib.__version__)\')"
    END_VERSIONS
    """
}
