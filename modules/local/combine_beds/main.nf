process COMBINE_BEDS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_upsetplot:0fc26c37f7821606' :
        'community.wave.seqera.io/library/polars_upsetplot:3382b69d3c1f6bf1' }"

    input:
    tuple val(meta), path(beds)
    val(max_shift)
    val(min_tools)
    val(min_samples)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: combined
    path "*.png"                                , emit: plots, optional: true
    path "*.json"                               , emit: multiqc, optional: true
    path "versions.yml"                         , emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    suffix      = task.ext.suffix      ?: "bed"
    template "combine.py"
}
