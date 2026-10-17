process COMBINEBEDS_FILTER {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_polars_pyarrow_upsetplot:8840b96e156438fc' :
        'community.wave.seqera.io/library/pandas_polars_pyarrow_upsetplot:6982d93f61d3e2ff' }"

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
    path "versions.yml"                         , emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    template "filter.py"

    stub:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        polars: \$(python -c "import polars; print(polars.__version__)")
        upsetplot: \$(python -c "import upsetplot; print(upsetplot.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
