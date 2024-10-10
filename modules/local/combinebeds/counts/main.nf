process COMBINEBEDS_COUNTS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_polars_pyarrow_upsetplot:8840b96e156438fc' :
        'community.wave.seqera.io/library/pandas_polars_pyarrow_upsetplot:6982d93f61d3e2ff' }"

    input:
    tuple val(meta), val(aggregation), path(candidates), path(beds)
    val(max_shift)
    val(consider_strand)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: combined, optional: true
    path "versions.yml"                         , emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    suffix      = task.ext.suffix      ?: "tsv"
    template "counts.py"
}
