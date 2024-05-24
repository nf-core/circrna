include { STAR_ALIGN as MATE1_1ST_PASS } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as MATE1_2ND_PASS } from '../../../modules/nf-core/star/align'
include { SJDB       as MATE1_SJDB     } from '../../../modules/local/star/sjdb'
include { STAR_ALIGN as MATE2_1ST_PASS } from '../../../modules/nf-core/star/align'
include { STAR_ALIGN as MATE2_2ND_PASS } from '../../../modules/nf-core/star/align'
include { SJDB       as MATE2_SJDB     } from '../../../modules/local/star/sjdb'
include { DCC        as MAIN           } from '../../../modules/local/dcc'
include { GAWK       as UNIFY          } from '../../../modules/nf-core/gawk'

workflow DCC {
    take:
    reads
    ch_fasta
    ch_gtf
    star_index
    star_junction
    ignore_sjdbgtf
    seq_platform
    seq_center
    bsj_reads
    
    main:
    ch_versions = Channel.empty()
    
    mate1 = reads.filter{ meta, reads -> !meta.single_end }
        .map{ meta, reads -> return [ [id: meta.id, single_end: true], reads[0] ] }
    MATE1_1ST_PASS( mate1, star_index, ch_gtf, ignore_sjdbgtf, seq_platform, seq_center )
    MATE1_SJDB( MATE1_1ST_PASS.out.tab
        .map{ meta, tab -> return tab }.collect().map{[[id: "mate1_sjdb"], it]}, bsj_reads )
    MATE1_2ND_PASS( mate1, star_index, MATE1_SJDB.out.sjtab, ignore_sjdbgtf, seq_platform, seq_center )

    mate2 = reads.filter{ meta, reads -> !meta.single_end }
        .map{ meta, reads -> return [ [id: meta.id, single_end: true], reads[1] ] }
    MATE2_1ST_PASS( mate2, star_index, ch_gtf, ignore_sjdbgtf, seq_platform, seq_center )
    MATE2_SJDB( MATE2_1ST_PASS.out.tab
        .map{ meta, tab -> return tab }.collect().map{[[id: "mate2_sjdb"], it]}, bsj_reads )
    MATE2_2ND_PASS( mate2, star_index, MATE2_SJDB.out.sjtab, ignore_sjdbgtf, seq_platform, seq_center )

    dcc_stage = star_junction.map{ meta, junction -> return [ meta.id, meta, junction]}
        .join(
            MATE1_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, junction] },
            remainder: true
        )
        .join(
            MATE2_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, junction] },
            remainder: true
        )
        .map{ id, meta, junction, mate1, mate2 -> return [ meta, junction, mate1, mate2 ]}

    dcc = dcc_stage.map{ it ->  [ it[0], it[1], it[2] ?: [], it[3] ?: [] ] }
    MAIN( dcc, ch_fasta.map{ meta, fasta -> fasta }, ch_gtf.map{ meta, gtf -> gtf } )
    UNIFY( MAIN.out.txt.map{ meta, txt -> [ meta + [tool: "dcc"], txt ] }, [] )

    ch_versions = ch_versions.mix(MATE1_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(MATE1_SJDB.out.versions)
    ch_versions = ch_versions.mix(MATE1_2ND_PASS.out.versions)
    ch_versions = ch_versions.mix(MATE2_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(MATE2_SJDB.out.versions)
    ch_versions = ch_versions.mix(MATE2_2ND_PASS.out.versions)
    ch_versions = ch_versions.mix(MAIN.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output
    
    versions = ch_versions
}