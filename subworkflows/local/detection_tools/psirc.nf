include { PSIRC_BSJ as BSJ } from '../../../modules/local/psirc/bsj'
include { PSIRC_FLI as FLI } from '../../../modules/local/psirc/fli'
include { GAWK as UNIFY    } from '../../../modules/nf-core/gawk'

workflow PSIRC {
    take:
    ch_reads
    ch_index
    detect_fli

    main:
    ch_versions = Channel.empty()

    BSJ(ch_reads, ch_index)
    ch_versions = ch_versions.mix(BSJ.out.versions)

    UNIFY(BSJ.out.bed, [], false)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    if (detect_fli) {
        FLI(
            ch_reads.join(BSJ.out.output),
            ch_index.map{ meta, transcriptome, _index -> [meta, transcriptome] }
        )
        ch_versions = ch_versions.mix(FLI.out.versions)
    }

    emit:
    bed = UNIFY.out.output.map{ meta, bed -> [meta + [tool: "psirc"], bed] }
    versions = ch_versions
}