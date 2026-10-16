process CIRIFULL_RO2 {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ciri-full:2.1.2--c949e1e8f8689713' :
        'community.wave.seqera.io/library/ciri-full:2.1.2--a656fc79dda2140f' }"

    input:
    tuple val(meta), path(sam), val(length)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}_ro2_info.list"), emit: list
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "2.1.2"
    """
    CIRI-full RO2 -r ${fasta} -s ${sam} -l ${length} -o ${prefix} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cirifull: ${VERSION}
    END_VERSIONS
    """
}
