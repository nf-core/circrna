process CIRCTOOLS_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/51/51f47a2e827b3bc883ac39107f50e25e5ad8f014138300106751c6a444dbff3d/data'
        : 'community.wave.seqera.io/library/gtfparse_polars:58e453d2bfabb473'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('gtfparse'), eval('python -c "import gtfparse; print(gtfparse.__version__)"'), emit: versions_gtfparse, topic: versions
    tuple val("${task.process}"), val('polars'), eval('python -c "import polars; print(polars.__version__)"'), emit: versions_polars, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('annotation.sh')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    """
}
