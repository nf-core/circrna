process CIRI3_UC_LIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/31cbdd6f827c754d8ebddcc66a9e50f97341b96b13c766c27cfd9c2b2f99f850/data'
:         'community.wave.seqera.io/library/python_polars_pyyaml:3d5318c405239fe0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${prefix}.ciri3.list"), emit: list
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'uc_list.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ciri3.list
    python - <<'PY'
    import platform
    import polars
    import yaml

    with open("versions.yml", "w") as versions:
        yaml.safe_dump({"${task.process}": {"python": platform.python_version(), "polars": polars.__version__, "pyyaml": yaml.__version__}}, versions)
    PY
    """
}
