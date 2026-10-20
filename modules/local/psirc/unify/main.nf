process PSIRC_UNIFY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b602845cfdd1b01ef0af6c42760522752eb1326c4eb74e8109277ad1156fe9bf/data'
        : 'community.wave.seqera.io/library/python_pandas:5febcd04d0f958b3'}"

    input:
    tuple val(meta), path(isoforms_tsv)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval('python --version | cut -d " " -f2'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas as pd; print(pd.__version__)"'), emit: versions_pandas, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'psirc_unify.sh'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    """
}
