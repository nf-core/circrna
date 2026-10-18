process CIRI3_CIRI3 {
    tag "${meta.id}"
    label 'process_high'

    conda 'environment.yml'
    container 'docker.io/cimg/openjdk@sha256:636ba84fbc06da072f6c48b0be7ef205de36eb1091f090c9b7255f6d95d8350c'

    input:
    tuple val(meta), path(sam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)
    path(circ_list)

    output:
    tuple val(meta), path("${prefix}.ciri3.txt"), emit: txt, optional: true
    tuple val(meta), path("${prefix}.ciri3.txt.BSJ_Matrix"), emit: bsj_matrix, optional: true
    tuple val(meta), path("${prefix}.ciri3.txt.FSJ_Matrix"), emit: fsj_matrix, optional: true
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def hasCircList = circ_list ? true : false
    def circListArg = hasCircList ? "-UC ${circ_list}" : ''
    // Limit the JVM heap to 80% of the task allocation to leave room for native memory.
    def heapMb = Math.max(1, (task.memory.toMega() * 0.8) as long)
    VERSION = '3.0.1'
    """
    java -Xmx${heapMb}m -jar \$(command -v CIRI3_Java_1.8.0.jar) \\
        -I ${sam} -O ${prefix}.ciri3.txt -F ${fasta} -A ${gtf} -T ${task.cpus} ${circListArg} ${args}

    if [ ${hasCircList} = false ] && [ ! -f ${prefix}.ciri3.txt ]; then
        echo 'CIRI3 detection did not create its required result file' >&2
        exit 1
    fi
    if [ ${hasCircList} = true ] && { [ ! -f ${prefix}.ciri3.txt.BSJ_Matrix ] || [ ! -f ${prefix}.ciri3.txt.FSJ_Matrix ]; }; then
        echo 'CIRI3 -UC mode did not create both required BSJ and FSJ matrices' >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri3: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = '3.0.1'
    """
    touch ${prefix}.ciri3.txt
    touch ${prefix}.ciri3.txt.BSJ_Matrix ${prefix}.ciri3.txt.FSJ_Matrix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri3: ${VERSION}
    END_VERSIONS
    """
}
