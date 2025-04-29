include { PSIRC_BSJ } from '../../../modules/local/psirc/bsj'

workflow PSIRC {
    take:
    ch_reads
    ch_index

    main:
    ch_versions = Channel.empty()

    PSIRC_BSJ(ch_reads, ch_index)
    ch_versions = ch_versions.mix(PSIRC_BSJ.out.versions)

    emit:
    bed = PSIRC_BSJ.out.bed
    versions = ch_versions
}