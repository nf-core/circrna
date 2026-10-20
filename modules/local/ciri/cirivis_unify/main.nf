process CIRIVIS_UNIFY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e43e750ee9bd4c15aab75f470900831786e8e4cd675708b3153373677b7882aa/data'
        : 'community.wave.seqera.io/library/pandas:3.0.6--73d876a9751f10fd'}"

    input:
    tuple val(meta), path(list)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('cirivis_unify.sh')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
