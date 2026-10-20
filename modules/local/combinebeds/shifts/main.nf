process COMBINEBEDS_SHIFTS {
    tag "${meta.id}"
    label "process_low"
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e241c2eab59c2ecaa0d7f44ad1441798cf129cb1ac178e492dd970c7fab8decb/data'
        : 'community.wave.seqera.io/library/python_altair_polars_vl-convert-python:bdb8533c65c6fe0b'}"

    input:
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("*.png"), emit: plots, optional: true
    tuple val(meta), path("*.json"), emit: multiqc, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: versions_polars, topic: versions
    tuple val("${task.process}"), val('altair'), eval("python -c 'import altair; print(altair.__version__)'"), emit: versions_altair, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template("shifts.sh")

    stub:
    """
    """
}
