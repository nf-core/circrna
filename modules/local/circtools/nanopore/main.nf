process CIRCTOOLS_NANOPORE {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/162d301bf9771dcc9bd647e7ea00cc9a959a933fc4b2bc0232fcb195acb09477/data'
        : 'community.wave.seqera.io/library/bedtools_circtools_nanofilt_pblat_samtools:20d0bb8a94c47e28'}"

    input:
    tuple val(meta), path(reads)
    path(reference)
    val(genome)

    output:
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p temp
    export TMPDIR="\$(pwd)/temp"

    circtools nanopore \
        -r \
        -s ${reads} \
        -R ${reference} \
        -C ${genome} \
        -O ${prefix} \
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circtools: \$(circtools -V)
    END_VERSIONS
    """
}
