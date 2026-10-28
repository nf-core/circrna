include { PSIRC_BSJ as BSJ         } from '../../../../modules/local/psirc/bsj'
include { PSIRC_BSJ_UNIFY as UNIFY } from '../../../../modules/local/psirc/bsj_unify'
include { INFER_STRAND             } from '../../../../modules/local/infer_strand'

workflow PSIRC {
    take:
    ch_reads
    ch_index
    ch_fasta
    ch_faidx

    main:

    BSJ(ch_reads, ch_index.first())

    // These are singleton references. Convert them to value channels so every sample task receives them.
    UNIFY(
        BSJ.out.bed.map { meta, bed -> [meta + [tool: "psirc"], bed] },
        ch_index.map { meta, transcriptome, _index -> [meta, transcriptome] }.first(),
    )

    INFER_STRAND(UNIFY.out.bed, ch_fasta.join(ch_faidx).first())

    emit:
    bed    = INFER_STRAND.out.bed
    output = BSJ.out.output
}
