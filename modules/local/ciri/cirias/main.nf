process CIRI_CIRIAS {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/ciri-full:2.1.2--c949e1e8f8689713'
        : 'community.wave.seqera.io/library/ciri-full:2.1.2--a656fc79dda2140f'}"

    input:
    tuple val(meta), path(sam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.0.0"
    """
    CIRI-AS -S ${sam} -C ${prefix}.ciri -F ${fasta} -A ${gtf} -O ${prefix} -D yes ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri-as: ${VERSION}
    END_VERSIONS
    """
}
