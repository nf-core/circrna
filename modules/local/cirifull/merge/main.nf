process CIRIFULL_MERGE {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/ciri-full:2.1.2--c949e1e8f8689713'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--a656fc79dda2140f'}"

    input:
    tuple val(meta), path(ciri), path(ciri_as), path(ro2)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}_merge_circRNA_detail.anno"), emit: anno
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "2.1.2"
    """
    CIRI-full Merge -a ${gtf} -r ${fasta} -c ${ciri} -as ${ciri_as} -ro ${ro2} -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cirifull: ${VERSION}
    END_VERSIONS
    """
}
