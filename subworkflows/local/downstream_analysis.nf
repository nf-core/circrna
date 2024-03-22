include { PREPARE_CIRCTEST } from '../../modules/local/circtest/prepare/main'

workflow DOWNSTREAM_ANALYSIS {
    take:
    linear_counts
    circ_counts

    main:
    ch_versions = Channel.empty()

    PREPARE_CIRCTEST(linear_counts, circ_counts)

    ch_versions = ch_versions.mix(PREPARE_CIRCTEST.out.versions)

    emit:
    versions = ch_versions
}