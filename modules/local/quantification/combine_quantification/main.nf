process COMBINE_QUANTIFICATION {
    label "process_single"

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
    tuple val(meta), path(inputs)
    tuple val(meta2), path(tx2gene)
    tuple val(meta3), path(circ_annotation)

    output:
    path("${meta.id}.linear.tsv")   , emit: linear
    path("${meta.id}.circular.tsv") , emit: circular
    path("${meta.id}.combined.tsv") , emit: combined
    path "versions.yml"             , emit: versions

    script:
    """
    combine_quantification.py --inputs ${inputs} \\
                                --tx2gene ${tx2gene} \\
                                --circ_annotation ${circ_annotation} \\
                                --out_linear ${meta.id}.linear.tsv \\
                                --out_circular ${meta.id}.circular.tsv \\
                                --out_combined ${meta.id}.combined.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
