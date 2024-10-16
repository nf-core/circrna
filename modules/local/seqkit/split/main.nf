process SEQKIT_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0' :
        'biocontainers/seqkit:2.8.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/*"), emit: split
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        split \\
        $args \\
        --threads $task.cpus \\
        $fasta \\
        --out-dir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
