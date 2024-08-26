include { CIRI_FULL } from './detection_tools/ciri_full'

workflow ISOFORM_DETECTION {

    take:
    ch_bwa_index
    ch_fasta
    ch_reads

    main:

    ch_versions = Channel.empty()

    CIRI_FULL ( ch_bwa_index, ch_fasta, ch_reads )
    ch_versions = ch_versions.mix(CIRI_FULL.out.versions)

    emit:

    versions = ch_versions
}
