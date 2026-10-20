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

    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('psirc_unify.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python --version 2>&1 | cut -d " " -f 2)"
        pandas: "\$(python -c \'import pandas; print(pandas.__version__)\')"
    END_VERSIONS
    """
}
