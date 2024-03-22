process PREPARE_CIRCTEST {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(gene_matrix)
    tuple val(meta2), path(circrna_matrix)

    output:
    tuple val(meta),  path("${meta.id}.circtest.tsv"), emit: genes
    tuple val(meta2), path("${meta2.id}.circtest.tsv"), emit: circs
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_circtest.py --in_genes ${gene_matrix} \\
                        --out_genes ${meta.id}.circtest.tsv \\
                        --in_circs ${circrna_matrix} \\
                        --out_circs ${meta2.id}.circtest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
