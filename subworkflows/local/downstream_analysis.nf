include { PREPARE_CIRCTEST } from '../../modules/local/circtest/prepare/main'
include { CIRCTEST         } from '../../modules/local/circtest/test/main'

workflow DOWNSTREAM_ANALYSIS {
    take:
    gene_counts
    circ_counts
    ch_phenotype

    main:
    ch_versions = Channel.empty()

    PREPARE_CIRCTEST(gene_counts, circ_counts)
    CIRCTEST(PREPARE_CIRCTEST.out.genes, PREPARE_CIRCTEST.out.circs, ch_phenotype)

    ch_versions = ch_versions.mix(PREPARE_CIRCTEST.out.versions)
    ch_versions = ch_versions.mix(CIRCTEST.out.versions)

    emit:
    versions = ch_versions
}