process CIRCTOOLS_RECONSTRUCT {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/circtools:2.0--b33a23b50d9f0697'
        : 'community.wave.seqera.io/library/circtools:2.0--f5bc60d7f93fefae'}"

    input:
    tuple val(meta), path(reads), path(bam)
    tuple val(meta2), path(annotation)

    output:
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    circtools reconstruct \
        --bamfile ${bam} \
        --annotation ${annotation} \
        --sampleName ${meta.id} \
        -C ${reads} \
        -T ./temp \
        -O ${prefix} \
        -P ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circtools: \$(circtools -V)
    END_VERSIONS
    """
}
