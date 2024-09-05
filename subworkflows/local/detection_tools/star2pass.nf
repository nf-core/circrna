include { STAR_ALIGN as PASS_1 } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as PASS_2 } from '../../../modules/nf-core/star/align'
include { SJDB                 } from '../../../modules/local/star/sjdb'


workflow STAR2PASS {
    take:
    reads
    star_index
    ch_gtf
    bsj_reads
    ignore_sjdbgtf
    seq_center
    seq_platform

    main:
    ch_versions = Channel.empty()

    PASS_1( reads, star_index, ch_gtf, ignore_sjdbgtf, seq_platform, seq_center)
    sjdb = PASS_1.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "star_sjdb"], it]}
    SJDB( sjdb, bsj_reads )
    PASS_2( reads, star_index, SJDB.out.sjtab, ignore_sjdbgtf, seq_platform, seq_center )

    ch_versions = ch_versions.mix(PASS_1.out.versions)
    ch_versions = ch_versions.mix(SJDB.out.versions)
    ch_versions = ch_versions.mix(PASS_2.out.versions)

    emit:
    junction = PASS_2.out.junction
    sam = PASS_2.out.sam
    tab = PASS_2.out.tab
    bam = PASS_2.out.bam

    versions = ch_versions
}
