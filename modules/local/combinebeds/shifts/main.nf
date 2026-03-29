process COMBINEBEDS_SHIFTS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb4ec5e7ad6feeda9b1cd6d194043b6cc0d2952c93a28cddc98fc7d67c078141/data' :
        'community.wave.seqera.io/library/altair_polars_pyyaml_vl-convert-python:c19053ed9a1a6146' }"

    input:
    tuple val(meta), path(beds)

    output:
    path "*.png"       , emit: plots, optional: true
    path "*.json"      , emit: multiqc, optional: true
    path "versions.yml", emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    template "shifts.py"
}
