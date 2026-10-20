process COMBINEBEDS_READS {
    tag "$meta.id"
    label "process_low"
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79a97132a439345fb983e26a8f44ffe1837a594eda02820d050b61bc02b8a96f/data'
        : 'community.wave.seqera.io/library/python_polars_seaborn_upsetplot_pruned:9bdb89c6afa3ce51'}"

    input:
    tuple val(meta), val(tools), path(beds)
    val(consider_strand)
    val(min_tools)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: combined
    path "*.png"                          , emit: plots, optional: true
    path "*.json"                         , emit: multiqc, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: versions_polars, topic: versions
    tuple val("${task.process}"), val('seaborn'), eval("python -c 'import seaborn; print(seaborn.__version__)'"), emit: versions_seaborn, topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval("python -c 'import matplotlib; print(matplotlib.__version__)'"), emit: versions_matplotlib, topic: versions
    tuple val("${task.process}"), val('upsetplot'), eval("python -c 'import upsetplot; print(upsetplot.__version__)'"), emit: versions_upsetplot, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template "reads.sh"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    """
}
