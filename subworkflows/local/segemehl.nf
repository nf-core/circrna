
include { SEGEMEHL_ALIGN } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_PARSE } from '../../modules/local/segemehl/parse/main'

workflow SEGEMEHL {

    take:
    reads
    fasta
    index
    bsj_reads

    main:
    ch_versions = Channel.empty()

    SEGEMEHL_ALIGN(
        reads,
        fasta,
        index
    )

    SEGEMEHL_PARSE(
        SEGEMEHL_ALIGN.out.results,
        bsj_reads
    )

    // Collect versions
    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)

    emit:
    segemehl_results = SEGEMEHL_PARSE.out.results.map{ meta, results -> meta.tool = "segemehl"; return [ meta, results ] }
    versions = ch_versions
}
