process COMBINE_BEDS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars:1.8.2--7ca02e3d227f4db6' :
        'community.wave.seqera.io/library/polars:1.8.2--2758744047b13d78' }"

    input:
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: combined
    path "versions.yml"                         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "bed"
    max_shift = task.ext.max_shift ?: 1
    min_files = task.ext.min_files ?: 1
    template "combine.py"
}
