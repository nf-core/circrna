process CIRCMINER_ORDERGTF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/49df43a275ec8361a3a4f810876e9a7e7034f9998084dccd8f77f299c18225cb/data'
:         'community.wave.seqera.io/library/python_gffutils_pyyaml:7521a2b68d4ddc13' }"

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

    python - <<'PY'
    import platform
    import gffutils
    import yaml

    with open("versions.yml", "w") as versions:
        yaml.safe_dump({"${task.process}": {"python": platform.python_version(), "gffutils": gffutils.__version__, "pyyaml": yaml.__version__}}, versions)
    PY
    """
}
