include { PSIRC_FLI as FLI     } from '../../../modules/local/psirc/fli'
include { PSIRC_UNIFY as UNIFY } from '../../../modules/local/psirc/unify'

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

    UNIFY(FLI.out.tsv)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    fasta = FLI.out.fasta.map{ meta, fasta -> [meta + [fli_tool: 'psirc'], fasta] }
    bed12 = UNIFY.out.bed.map{ meta, bed -> [meta + [fli_tool: 'psirc'], bed] }

    versions = ch_versions
}
