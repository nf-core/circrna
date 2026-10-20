process MAJORITY_VOTE {
    tag "$meta.id"
    label 'process_medium'
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c0e5925c40ad9479f2280efd1f3719c4ca7e1afe35291e8aa77fa47f79ce0d1c/data'
        : 'community.wave.seqera.io/library/majority_vote:4af2b9d45db454df'}"

    input:
    tuple val(meta), path(bindingsites)
    val(min_tools)

    output:
    tuple val(meta), path("${meta.id}.majority.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.targets.tsv") , emit: targets
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: versions_polars, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'majority.sh'

    stub:
    """
    touch ${meta.id}.majority.tsv
    touch ${meta.id}.targets.tsv

    """
}
