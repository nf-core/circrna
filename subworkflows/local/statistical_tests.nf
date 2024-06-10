include { CIRCTEST_PREPARE } from '../../modules/local/circtest/prepare'

workflow STATISTICAL_TESTS {
    take:
    ch_quantification
    ch_gene_counts
    ch_circ_counts
    ch_phenotype
    
    main:
    ch_versions = Channel.empty()

    CIRCTEST_PREPARE(ch_circ_counts, ch_gene_counts)
}