process COMBINE_BEDS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_upsetplot:0fc26c37f7821606' :
        'community.wave.seqera.io/library/polars_upsetplot:3382b69d3c1f6bf1' }"

    input:
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: combined
    path "*.png"                                , emit: plots
    path "*.json"                               , emit: multiqc
    path "versions.yml"                         , emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    suffix      = task.ext.suffix      ?: "bed"
    max_shift   = task.ext.max_shift   ?: 1
    min_tools   = task.ext.min_tools   ?: 1
    min_samples = task.ext.min_samples ?: 1
    template "combine.py"
}
