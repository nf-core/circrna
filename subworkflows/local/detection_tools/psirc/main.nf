include { PSIRC_BSJ as BSJ } from '../../../../modules/local/psirc/bsj'
include { PSIRC_BSJ_UNIFY as UNIFY } from '../../../../modules/local/psirc/bsj_unify'
include { INFER_STRAND } from '../../../../modules/local/infer_strand'

workflow PSIRC {
    take:
    ch_reads
    ch_index
    ch_fasta
    ch_faidx

    main:
    ch_versions = channel.empty()

    BSJ(ch_reads, ch_index)
    ch_versions = ch_versions.mix(BSJ.out.versions)

    // These are singleton references. Convert them to value channels so every sample task receives them.
    UNIFY(
        BSJ.out.bed.map { meta, bed -> [meta + [tool: "psirc"], bed] },
        ch_index.map { meta, transcriptome, _index -> [meta, transcriptome] }.first(),
    )
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    INFER_STRAND(UNIFY.out.bed, ch_fasta.join(ch_faidx).first())
    ch_versions = ch_versions.mix(INFER_STRAND.out.versions)

    emit:
    bed      = INFER_STRAND.out.bed
    output   = BSJ.out.output

    versions = ch_versions
}
