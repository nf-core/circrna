process COMBINEBEDS_READS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_pyyaml_seaborn_upsetplot:f266bc0cdfb4cc77' :
        'community.wave.seqera.io/library/polars_pyyaml_seaborn_upsetplot:9024174fe22d5591' }"

    input:
    tuple val(meta), val(tools), path(beds)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: combined
    path "*.png"                          , emit: plots, optional: true
    path "*.json"                         , emit: multiqc, optional: true
    path "versions.yml"                   , emit: versions

    script:
    prefix          = task.ext.prefix          ?: "${meta.id}"
    consider_strand = task.ext.consider_strand ?: false
    min_tools       = task.ext.min_tools       ?: 2
    template "reads.py"
}
