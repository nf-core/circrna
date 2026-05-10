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

    main:
    ch_versions = channel.empty()

    PASS_1( reads, star_index, ch_gtf, ignore_sjdbgtf)
    sjdb = PASS_1.out.tab.map{ _meta, tab -> return tab }.collect().map{ tabs -> [[id: "star_sjdb"], tabs]}

    SJDB( sjdb, bsj_reads )
    ch_versions = ch_versions.mix(SJDB.out.versions)

    PASS_2( reads, star_index, SJDB.out.sjtab, ignore_sjdbgtf )

    emit:
    junction = PASS_2.out.junction
    sam = PASS_2.out.sam
    tab = PASS_2.out.tab
    bam = PASS_2.out.bam

    versions = ch_versions
}
