include { CIRCMINER_CIRCMINER as DETECT } from '../../../../modules/local/circminer/circminer'
include { GAWK as UNIFY                 } from '../../../../modules/nf-core/gawk'
include { INFER_STRAND                  } from '../../../../modules/local/infer_strand'

workflow CIRCMINER {
    take:
    reads
    index
    gtf
    fasta
    faidx

    main:

    DETECT(reads, index.first(), gtf.first())

    UNIFY(
        DETECT.out.report.map { meta, report -> [meta + [tool: 'circminer'], report, 'bed'] },
        [],
        false,
    )
    ch_reference = fasta
        .join(faidx)
        .first()

    INFER_STRAND(UNIFY.out.output, ch_reference)

    emit:
    bed = INFER_STRAND.out.bed
}
