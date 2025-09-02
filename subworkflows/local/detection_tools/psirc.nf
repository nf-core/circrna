include { PSIRC_BSJ as BSJ } from '../../../modules/local/psirc/bsj'
include { GAWK as UNIFY    } from '../../../modules/nf-core/gawk'

workflow PSIRC {
    take:
    ch_reads
    ch_index

    main:
    ch_versions = Channel.empty()

    BSJ(ch_reads, ch_index)
    ch_versions = ch_versions.mix(BSJ.out.versions)

    UNIFY(
        BSJ.out.bed.map { meta, bed -> [meta + [tool: "psirc"], bed] },
        [],
        false,
    )
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed      = UNIFY.out.output
    output   = BSJ.out.output

    versions = ch_versions
}
