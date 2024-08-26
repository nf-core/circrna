include { CIRI_FULL } from './detection_tools/ciri_full'

workflow ISOFORM_DETECTION {

    take:
    ch_reads
    ch_reads_untrimmed
    ch_fasta
    ch_gtf

    main:

    ch_versions = Channel.empty()

    CIRI_FULL ( ch_reads_untrimmed, ch_fasta, ch_gtf )
    ch_versions = ch_versions.mix(CIRI_FULL.out.versions)

    emit:

    versions = ch_versions
}
