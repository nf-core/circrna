process CIRI3_UC_MATRIX {
    tag "$meta2.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/31cbdd6f827c754d8ebddcc66a9e50f97341b96b13c766c27cfd9c2b2f99f850/data'
:         'community.wave.seqera.io/library/python_polars_pyyaml:3d5318c405239fe0' }"

    input:
    tuple val(meta), path(circ_list)
    tuple val(meta2), path(bsj_matrix), path(fsj_matrix)

    output:
    tuple val(meta2), path("${prefix}.bsj.tsv"), emit: bsj
    tuple val(meta2), path("${prefix}.fsj.tsv"), emit: fsj
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta2.id}"
    template 'uc_matrix.py'

    stub:
    prefix = task.ext.prefix ?: "${meta2.id}"
    """
    printf 'circRNA_ID\\tcount\\n' > ${prefix}.bsj.tsv
    printf 'circRNA_ID\\tcount\\n' > ${prefix}.fsj.tsv
    python - <<'PY'
    import platform
    import polars
    import yaml

    with open("versions.yml", "w") as versions:
        yaml.safe_dump({"${task.process}": {"python": platform.python_version(), "polars": polars.__version__, "pyyaml": yaml.__version__}}, versions)
    PY
    """
}
