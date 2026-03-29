include { CIRCTOOLS_LONGREAD } from './longread_tools/circtools'

workflow LONGREAD {
    take:
    ch_reads
    genome

    main:
    ch_versions = Channel.empty()

    CIRCTOOLS_LONGREAD(ch_reads, genome)
    ch_versions = ch_versions.mix(CIRCTOOLS_LONGREAD.out.versions)

    emit:
    versions = ch_versions
}
