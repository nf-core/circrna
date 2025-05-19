include { PSIRC_FLI as FLI } from '../../../modules/local/psirc/fli'

workflow PSIRC {
    take:
    ch_reads
    ch_psirc_bsj
    ch_psirc_index

    main:
    ch_versions = Channel.empty()

    FLI(
        ch_reads.join(ch_psirc_bsj),
        ch_psirc_index.map { meta, transcriptome, _index -> [meta, transcriptome] },
    )
    ch_versions = ch_versions.mix(FLI.out.versions)

    emit:
    versions = ch_versions
}
