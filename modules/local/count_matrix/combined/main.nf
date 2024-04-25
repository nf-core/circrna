process COUNTS_COMBINED {
    label "process_single"

    conda "bioconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(beds)

    output:
    path("merged_counts.bed") , emit: counts_bed
    path("merged_counts.tsv") , emit: counts_tsv
    path "versions.yml"       , emit: versions

    script:
    """
    counts_combined.py --beds ${beds} --out_bed merged_counts.bed --out_tsv merged_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
