include { PSIRC_BSJ } from '../../../modules/local/psirc/bsj'
include { PSIRC_FLI } from '../../../modules/local/psirc/fli'

workflow PSIRC {
    take:
    ch_reads
    ch_index

    main:
    ch_versions = Channel.empty()

    PSIRC_BSJ(ch_reads, ch_index)
    ch_versions = ch_versions.mix(PSIRC_BSJ.out.versions)

    PSIRC_FLI(
        ch_reads.join(PSIRC_BSJ.out.output),
        ch_index.map{ meta, transcriptome, _index -> [meta, transcriptome] }
    )
    ch_versions = ch_versions.mix(PSIRC_FLI.out.versions)

    emit:
    bed = PSIRC_BSJ.out.bed
    versions = ch_versions
}