include { SEGEMEHL_ALIGN              } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_PARSE              } from '../../modules/local/segemehl/parse/main'
include { STAR_ALIGN as STAR_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_2ND_PASS } from '../../modules/nf-core/star/align/main'

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

    segemehl_parse = SEGEMEHL_ALIGN.out.results.map{ meta, results ->  meta.tool = "segemehl"; return [ meta, results ] }

    SEGEMEHL_PARSE(
        segemehl_parse,
        bsj_reads
    )

    STAR_1ST_PASS(
        reads,
        star_index,
        gtf,
        true,
        '',
        ''
    )

    STAR_2ND_PASS(
        reads,
        star_index,
        STAR_1ST_PASS.out.junction.collect(), // use Chimeric Junctions in place of GTF, implement in modules.config.
        true,
        '',
        ''
    )

    // collect versions
    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)
    ch_versions = ch_versions.mix(STAR_1ST_PASS.out.versions)

    emit:
    segemehl_results = SEGEMEHL_PARSE.out.results
    versions = ch_versions
}
