process MAJORITY_VOTE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
            'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(bindingsites)


    output:
    tuple val(meta), path("${meta.id}.majority.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.targets.tsv") , emit: targets
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
