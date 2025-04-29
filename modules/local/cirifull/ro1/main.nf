process CIRIFULL_RO1 {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ciri-full:2.1.2--c949e1e8f8689713' :
        'community.wave.seqera.io/library/ciri-full:2.1.2--a656fc79dda2140f' }"

    input:
    tuple val(meta), path(reads)

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI-full RO1 -1 ${reads[0]} -2 ${reads[1]} -o ${prefix} $args
    """
}