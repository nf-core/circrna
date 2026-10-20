process COMBINEBEDS_FILTER {
    tag "$meta.id"
    label "process_low"
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/48af262bc66e75295b1a314cbf45b01bded343911ffeaa0357009f7370cb4573/data'
        : 'community.wave.seqera.io/library/python_polars_upsetplot_pandas_matplotlib:a16ba3398a0ae7a3'}"

    input:
    tuple val(meta), path(beds)
    val(max_shift)
    val(consider_strand)
    val(min_tools)
    val(min_samples)
    val(suffix)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: combined, optional: true
    path "*.png"                                , emit: plots, optional: true
    path "*.json"                               , emit: multiqc, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: versions_polars, topic: versions
    tuple val("${task.process}"), val('upsetplot'), eval("python -c 'import upsetplot; print(upsetplot.__version__)'"), emit: versions_upsetplot, topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval("python -c 'import matplotlib; print(matplotlib.__version__)'"), emit: versions_matplotlib, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template "filter.sh"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${suffix}

    """
}
