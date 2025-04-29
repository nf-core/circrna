process CIRI_CIRI2 {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ciri-full:2.1.2--c949e1e8f8689713' :
        'community.wave.seqera.io/library/ciri-full:2.1.2--a656fc79dda2140f' }"

    input:
    tuple val(meta), path(sam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI -I $sam -O $prefix -F $fasta -A $gtf -T $task.cpus
    """
}