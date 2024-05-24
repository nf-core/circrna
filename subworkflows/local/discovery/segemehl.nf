include { SEGEMEHL_INDEX as INDEX   } from '../../../modules/nf-core/segemehl/index'
include { SEGEMEHL_ALIGN as ALIGN   } from '../../../modules/nf-core/segemehl/align'
include { SEGEMEHL_FILTER as FILTER } from '../../../modules/local/segemehl/filter'


workflow SEGEMEHL {
    take:
    reads
    fasta
    index
    bsj_reads

    main:
    ch_versions = Channel.empty()

    index = index ?: INDEX( fasta ).index
    ALIGN( reads, fasta, index )
    FILTER( ALIGN.out.results
        .map{ meta, results ->  [ meta + [tool: "segemehl"], results ] }, bsj_reads )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(FILTER.out.versions)

    emit:
    matrix = FILTER.out.matrix
    results = FILTER.out.results

    versions = ch_versions
}