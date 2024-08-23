include { CIRCTEST_PREPARE  } from '../../modules/local/circtest/prepare'
include { CIRCTEST_CIRCTEST } from '../../modules/local/circtest/circtest'

workflow STATISTICAL_TESTS {
    take:
    ch_gene_counts
    ch_circ_counts
    ch_phenotype

    main:
    ch_versions = Channel.empty()

    ch_counts = ch_gene_counts.join(ch_circ_counts)

    CIRCTEST_PREPARE(ch_counts)
    ch_versions = ch_versions.mix(CIRCTEST_PREPARE.out.versions)

    CIRCTEST_CIRCTEST(CIRCTEST_PREPARE.out.counts,
        ch_phenotype)
    ch_versions = ch_versions.mix(CIRCTEST_CIRCTEST.out.versions)

    emit:
    versions = ch_versions
}
