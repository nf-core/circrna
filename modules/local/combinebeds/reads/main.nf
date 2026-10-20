process COMBINEBEDS_READS {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79a97132a439345fb983e26a8f44ffe1837a594eda02820d050b61bc02b8a96f/data'
        : 'community.wave.seqera.io/library/python_polars_seaborn_upsetplot_pruned:9bdb89c6afa3ce51'}"

    input:
    tuple val(meta), val(tools), path(beds)
    val consider_strand
    val min_tools

    output:
    tuple val(meta), path("${prefix}.bed"), emit: combined
    path "*.png", emit: plots, optional: true
    path "*.json", emit: multiqc, optional: true

    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template("reads.py")

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python --version 2>&1 | cut -d " " -f 2)"
        polars: "\$(python -c \'import polars; print(polars.__version__)\')"
        seaborn: "\$(python -c \'import seaborn; print(seaborn.__version__)\')"
        matplotlib: "\$(python -c \'import matplotlib; print(matplotlib.__version__)\')"
        upsetplot: "\$(python -c \'import upsetplot; print(upsetplot.__version__)\')"
    END_VERSIONS
    """
}
