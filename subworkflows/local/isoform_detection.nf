include { CIRI_FULL } from './isoform_tools/ciri_full'

workflow ISOFORM_DETECTION {

    take:
    ch_reads
    ch_reads_untrimmed
    ch_fasta
    ch_gtf

    main:

    ch_versions = Channel.empty()

    tools_selected = params.tools.split(',').collect{it.trim().toLowerCase()}

    if (tools_selected.contains('ciri-full')) {
        CIRI_FULL ( ch_reads_untrimmed, ch_fasta, ch_gtf )
        ch_versions = ch_versions.mix(CIRI_FULL.out.versions)
    }

    emit:

    versions = ch_versions
}
