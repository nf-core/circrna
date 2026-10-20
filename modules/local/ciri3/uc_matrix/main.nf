process CIRI3_UC_MATRIX {
    tag "${meta2.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc24a035d2b10ae5cddc290f972da775bed663c0a9f3d046386412b2922112ad/data'
        : 'community.wave.seqera.io/library/python_polars:708bbfa3590a8658'}"

    input:
    tuple val(meta), path(circ_list)
    tuple val(meta2), path(bsj_matrix), path(fsj_matrix)

    output:
    tuple val(meta2), path("${prefix}.bsj.tsv"), emit: bsj
    tuple val(meta2), path("${prefix}.fsj.tsv"), emit: fsj
    tuple val("${task.process}"), val('python'), eval('python --version | cut -d " " -f2'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval('python -c "import polars; print(polars.__version__)"'), emit: versions_polars, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta2.id}"
    template('uc_matrix.sh')

    stub:
    prefix = task.ext.prefix ?: "${meta2.id}"
    """
    printf 'circRNA_ID\\tcount\\n' > ${prefix}.bsj.tsv
    printf 'circRNA_ID\\tcount\\n' > ${prefix}.fsj.tsv
    """
}
