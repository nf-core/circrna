include { SEGEMEHL_INDEX as INDEX   } from '../../../modules/nf-core/segemehl/index'
include { SEGEMEHL_ALIGN as ALIGN   } from '../../../modules/nf-core/segemehl/align'
include { GAWK as EXTRACT           } from '../../../modules/nf-core/gawk'
include { GNU_SORT as SORT          } from '../../../modules/nf-core/gnu/sort'
include { BEDTOOLS_GROUPBY as GROUP } from '../../../modules/nf-core/bedtools/groupby'
include { GAWK as UNIFY             } from '../../../modules/nf-core/gawk'

workflow SEGEMEHL {
    take:
    reads
    fasta
    index

    main:
    ch_versions = Channel.empty()

    index = index ?: INDEX( fasta ).index

    ALIGN( reads, fasta, index )
    EXTRACT( ALIGN.out.single_bed
        .map{ meta, bed ->  [ meta + [tool: "segemehl"], bed ] }, [] )

    SORT( EXTRACT.out.output )
    GROUP( SORT.out.sorted, 5 )
    UNIFY( GROUP.out.bed, [] )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(EXTRACT.out.versions)
    ch_versions = ch_versions.mix(SORT.out.versions)
    ch_versions = ch_versions.mix(GROUP.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}
