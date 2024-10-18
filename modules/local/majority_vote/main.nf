process MAJORITY_VOTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'oras://community.wave.seqera.io/library/polars_pyyaml:962a0cf7480258c7' :
            'community.wave.seqera.io/library/polars_pyyaml:ad93db0d7bcd508e' }"

    input:
    tuple val(meta), path(bindingsites)

    output:
    tuple val(meta), path("${meta.id}.majority.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.targets.tsv") , emit: targets
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    min_tools = params.mirna_min_tools
    template 'majority.py'

    stub:
    """
    touch ${meta.id}.majority.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
