process SPLIT_TYPES {
    tag "${meta.id}"
    label 'process_single'
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fadc2c277b689c5c075f8ca06347864b19cb30d5355f819cba851a3c24397f77/data'
        : 'community.wave.seqera.io/library/split_types:84ca943aca0dd262'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("linear.tsv"), emit: linear
    tuple val(meta), path("circular.tsv"), emit: circular
    tuple val("${task.process}"), val('gawk'), eval("gawk --version | sed -n '1s/,.*//p' | cut -d ' ' -f 3"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gawk -F'\\t' \\
        'NR==1 {print > "circular.tsv"; print > "linear.tsv"} \\
        NR>1 {if (\$1 ~ /^.+:[0-9]+-[0-9]+(:[+-])?/) print > "circular.tsv"; else print > "linear.tsv"}' ${input}
    """

    stub:
    """
    touch linear.tsv
    touch circular.tsv
    """
}
