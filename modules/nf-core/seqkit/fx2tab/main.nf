process SEQKIT_FX2TAB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c854c28b0f0a7a895f283f230367875922689c9b20fa1bca5474bc9999e56ff2/data'
        : 'community.wave.seqera.io/library/seqkit:2.14.0--fe9cf2addeefb114'}"

    input:
    tuple val(meta), path(fastx)
    val out_ext

    output:
    tuple val(meta), path("*.txt*"), emit: text
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = out_ext ?: "txt.zst"
    """
    seqkit \\
        fx2tab \\
        ${args} \\
        --threads ${task.cpus} \\
        ${fastx} \\
        -o ${prefix}.${suffix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = out_ext ?: "txt.zst"
    """
    echo ${args}

    touch ${prefix}.${suffix}
    """
}
