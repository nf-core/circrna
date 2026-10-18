process CIRCMINER_ORDERGTF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/18/1841daa69f98a0b0ffcb8f545070c8350a75febb167202136eab0990131d31c0/data'
:         'community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927' }"

    input:
    tuple val(meta), path(gtf, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
    path 'versions.yml'                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.circminer"
    template 'ordergtf.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.circminer"
    """
    touch ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
