include { CIRI3_CIRI3 as CIRI3_RUN } from '../../../../modules/local/ciri3/ciri3'
include { GAWK as UNIFY } from '../../../../modules/nf-core/gawk'

workflow CIRI3 {
    take:
    ch_sam
    ch_fasta
    ch_gtf

    main:
    ch_versions = channel.empty()

    CIRI3_RUN(ch_sam, ch_fasta.first(), ch_gtf.first(), [])
    ch_versions = ch_versions.mix(CIRI3_RUN.out.versions)

    UNIFY(
        CIRI3_RUN.out.txt.map { meta, txt -> [meta + [tool: 'ciri3'], txt, 'ciri3.bed'] },
        [],
        false,
    )

    emit:
    bed = UNIFY.out.output
    versions = ch_versions
}
