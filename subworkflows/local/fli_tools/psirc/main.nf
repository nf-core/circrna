include { PSIRC_FLI as FLI     } from '../../../../modules/local/psirc/fli'
include { PSIRC_UNIFY as UNIFY } from '../../../../modules/local/psirc/unify'

workflow PSIRC {
    take:
    ch_reads
    ch_psirc_bsj
    ch_psirc_index

    main:

    FLI(
        ch_reads.join(ch_psirc_bsj),
        ch_psirc_index.map { meta, transcriptome, _index -> [meta, transcriptome] },
    )

    UNIFY(FLI.out.tsv)

    emit:
    fasta = FLI.out.fasta.map { meta, fasta -> [meta + [fli_tool: 'psirc'], fasta] }
    bed12 = UNIFY.out.bed.map { meta, bed -> [meta + [fli_tool: 'psirc'], bed] }
}
