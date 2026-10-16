process CIRCTOOLS_REFERENCE {
    tag "${genome}"
    label 'process_low'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/circtools:2.0--b33a23b50d9f0697'
        : 'community.wave.seqera.io/library/circtools:2.0--f5bc60d7f93fefae'}"

    input:
    val(genome)

    output:
    path "reference", emit: reference
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p temp
    export TMPDIR=./temp
    circtools nanopore -d -R reference/ -C ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circtools: \$(circtools -V)
    END_VERSIONS
    """
}
