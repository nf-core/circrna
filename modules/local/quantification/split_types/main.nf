process SPLIT_TYPES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("linear.tsv")  , emit: linear
    tuple val(meta), path("circular.tsv"), emit: circular
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -F'\\t' \\
        'NR==1 {print > "circular.tsv"; print > "linear.tsv"} \\
        NR>1 {if (\$1 ~ /^circ_/) print > "circular.tsv"; else print > "linear.tsv"}' ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    """
    touch linear.tsv
    touch circular.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
