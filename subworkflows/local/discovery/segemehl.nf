include { SEGEMEHL_INDEX as INDEX   } from '../../../modules/nf-core/segemehl/index'
include { SEGEMEHL_ALIGN as ALIGN   } from '../../../modules/nf-core/segemehl/align'
include { GAWK as UNIFY             } from '../../../modules/nf-core/gawk'
include { BEDTOOLS_GROUPBY as GROUP } from '../../../modules/nf-core/bedtools/groupby'

workflow SEGEMEHL {
    take:
    reads
    fasta
    index

    main:
    ch_versions = Channel.empty()

    index = index ?: INDEX( fasta ).index
    ALIGN( reads, fasta, index )
    UNIFY( ALIGN.out.single_bed
        .map{ meta, bed ->  [ meta + [tool: "segemehl"], bed ] }, [] )

    GROUP( UNIFY.out.output, 5 )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)
    ch_versions = ch_versions.mix(GROUP.out.versions)

    emit:
    bed = GROUP.out.bed

    versions = ch_versions
}
