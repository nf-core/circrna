include { CIRCTOOLS_REFERENCE } from '../../../modules/local/circtools/reference'
include { CIRCTOOLS_NANOPORE } from '../../../modules/local/circtools/nanopore'

workflow CIRCTOOLS_LONGREAD {
    take:
    ch_reads
    genome

    main:
    ch_versions = Channel.empty()

    CIRCTOOLS_REFERENCE(genome)
    ch_versions = ch_versions.mix(CIRCTOOLS_REFERENCE.out.versions)

    CIRCTOOLS_NANOPORE(ch_reads, CIRCTOOLS_REFERENCE.out.reference, genome)
    ch_versions = ch_versions.mix(CIRCTOOLS_NANOPORE.out.versions)

    emit:
    versions = ch_versions
}
