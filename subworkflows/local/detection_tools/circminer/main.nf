include { CIRCMINER_CIRCMINER as DETECT } from '../../../../modules/local/circminer/circminer'
include { CIRCMINER_PREPAREGTF as PREPAREGTF } from '../../../../modules/local/circminer/preparegtf'
include { GAWK as UNIFY } from '../../../../modules/nf-core/gawk'
include { INFER_STRAND } from '../../../../modules/local/infer_strand'

workflow CIRCMINER {
    take:
    reads
    index
    gtf
    fasta
    faidx

    main:
    ch_versions = channel.empty()

    PREPAREGTF(gtf.first())
    ch_versions = ch_versions.mix(PREPAREGTF.out.versions)

    DETECT(reads, index.first(), PREPAREGTF.out.gtf.first())
    ch_versions = ch_versions.mix(DETECT.out.versions)

    UNIFY(
        DETECT.out.report.map { meta, report -> [meta + [tool: 'circminer'], report, 'bed'] },
        [],
        false,
    )
    ch_reference = fasta
        .join(faidx)
        .first()

    INFER_STRAND(UNIFY.out.output, ch_reference)
    ch_versions = ch_versions.mix(INFER_STRAND.out.versions)

    emit:
    bed = INFER_STRAND.out.bed
    versions = ch_versions
}
