process CIRCTEST_PREPARE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/r-base:4.2.1"

    input:
    tuple val(meta), path(gene_counts), path(circ_counts)

    output:
    tuple val(meta), path('*_genes.tsv'), path('*_circs.tsv'), emit: counts, optional: true

    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'prepare.R'
}
