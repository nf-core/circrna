process SEQKIT_SPLIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5266dd21f30e777f3b05d029e5b523e5a3f78ff2659cc23d58f2175700c20815/data'
        : 'community.wave.seqera.io/library/seqkit_split:b24e98aca17fcd40'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/*"), emit: split
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^seqkit v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        split \\
        ${args} \\
        --threads ${task.cpus} \\
        ${fasta} \\
        --out-dir ${prefix}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${fasta.baseName}.part_001.${fasta.extension}

    """
}
