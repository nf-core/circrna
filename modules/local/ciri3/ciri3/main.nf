process CIRI3_CIRI3 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a64fe42a183f0972414bddfa15d354a8da34b16e120ca984064d13d881aa2fbf/data'
        : 'community.wave.seqera.io/library/ciri3:3.0.1--9d1e9d05ff71b913'}"

    input:
    tuple val(meta), path(sam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)
    path circ_list

    output:
    tuple val(meta), path("${prefix}.ciri3.txt"), emit: txt, optional: true
    tuple val(meta), path("${prefix}.ciri3.txt.BSJ_Matrix"), emit: bsj_matrix, optional: true
    tuple val(meta), path("${prefix}.ciri3.txt.FSJ_Matrix"), emit: fsj_matrix, optional: true
    // CIRI3 has no version option; this is the pinned Bioconda package version.
    tuple val("${task.process}"), val('ciri3'), val('3.0.1'), emit: versions_ciri3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def hasCircList = circ_list ? true : false
    def circListArg = hasCircList ? "-UC ${circ_list}" : ''
    // Limit the JVM heap to 80% of the task allocation to leave room for native memory.
    def heapMb = Math.max(1, (task.memory.toMega() * 0.8) as long)
    """
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:+\$JAVA_TOOL_OPTIONS }-Xmx${heapMb}m"
    ciri3 \\
        -I ${sam} -O ${prefix}.ciri3.txt -F ${fasta} -A ${gtf} -T ${task.cpus} ${circListArg} ${args}

    if [ ${hasCircList} = false ] && [ ! -f ${prefix}.ciri3.txt ]; then
        echo 'CIRI3 detection did not create its required result file' >&2
        exit 1
    fi
    if [ ${hasCircList} = true ] && { [ ! -f ${prefix}.ciri3.txt.BSJ_Matrix ] || [ ! -f ${prefix}.ciri3.txt.FSJ_Matrix ]; }; then
        echo 'CIRI3 -UC mode did not create both required BSJ and FSJ matrices' >&2
        exit 1
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ciri3.txt
    touch ${prefix}.ciri3.txt.BSJ_Matrix ${prefix}.ciri3.txt.FSJ_Matrix
    """
}
