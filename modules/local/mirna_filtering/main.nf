process MIRNA_FILTERING {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(normalized_counts)
    val(mirna_min_sample_percentage)
    val(mirna_min_reads)

    output:
    tuple val(meta), path("${meta.id}.normalized_counts_filtered.tsv"), emit: filtered
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'mirna_filtering.R'

    stub:
    """
    touch ${meta.id}.normalized_counts_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1s/^R version \\([0-9.]*\\).*\$/\\1/p')
    END_VERSIONS
    """
}
