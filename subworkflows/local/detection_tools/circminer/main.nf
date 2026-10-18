include { CIRCMINER_CIRCMINER as DETECT } from '../../../../modules/local/circminer/circminer'
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

    DETECT(reads, index, gtf)
    ch_versions = ch_versions.mix(DETECT.out.versions)

    UNIFY(
        DETECT.out.report.map { meta, report -> [meta + [tool: 'circminer'], report, 'bed'] },
        [],
        false,
    )
    ch_reference = fasta
        .join(faidx)
        .map { meta, reference, reference_fai -> [meta, reference, reference_fai] }

    INFER_STRAND(UNIFY.out.output, ch_reference)
    ch_versions = ch_versions.mix(INFER_STRAND.out.versions)

    emit:
    bed = INFER_STRAND.out.bed
    versions = ch_versions
}
