include { SEGEMEHL_ALIGN } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_PARSE } from '../../modules/local/segemehl/parse/main'

workflow CIRCRNA_DISCOVERY {

    take:
    reads
    fasta
    gtf
    bowtie_index
    bowtie2_index
    segemehl_index
    star_index
    bsj_reads

    main:
    ch_versions = Channel.empty()

    SEGEMEHL_ALIGN(
        reads,
        fasta,
        segemehl_index
    )

    // append tool name to meta
    segemehl_parse = SEGEMEHL_ALIGN.out.results.map{ meta, results ->  meta.tool = "segemehl"; return [ meta, results ] }

    SEGEMEHL_PARSE(
        segemehl_parse,
        bsj_reads
    )

    // Collect versions
    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)

    emit:
    segemehl_results = SEGEMEHL_PARSE.out.results
    versions = ch_versions
}
