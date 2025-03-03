process COMBINEBEDS_SHIFTS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/altair_polars_vl-convert-python:e6f1dca28de76d13' :
        'community.wave.seqera.io/library/altair_polars_vl-convert-python:a6c5ee679445250d' }"

    input:
    tuple val(meta), path(beds)

    output:
    path "*.png"       , emit: plots
    path "*.json"      , emit: multiqc
    path "versions.yml", emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    template "shifts.py"
}
