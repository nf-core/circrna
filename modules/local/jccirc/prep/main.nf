process JCCIRC_PREP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc8e5ab4ba32d58b69d5395ae4b39f0fde2467e5a81c8a6895d664309bb27498/data'
        : 'community.wave.seqera.io/library/polars:2.0.0--2895a858bbcef9a5'}"

    input:
    tuple val(meta), path(bsj_annotation), path(bsj_reads)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: merged
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval('python -c "import polars; print(polars.__version__)"'), emit: versions_polars, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('prep.sh')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    """
}
