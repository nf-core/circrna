
include { ANNOTATION                  } from '../../modules/local/annotation/main'
include { CIRCEXPLORER2_REFERENCE     } from '../../modules/local/circexplorer2/reference/main'
include { CIRCEXPLORER2_PARSE         } from '../../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE      } from '../../modules/nf-core/circexplorer2/annotate/main'
include { CIRCEXPLORER2_FILTER        } from '../../modules/local/circexplorer2/filter/main'
include { CIRCRNA_FINDER_FILTER       } from '../../modules/local/circrna_finder/filter/main'
include { SEGEMEHL_ALIGN              } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_FILTER             } from '../../modules/local/segemehl/filter/main'
include { STAR_ALIGN as STAR_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB                        } from '../../modules/local/star/sjdb/main'

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

    //
    // SEGEMEHL WORKFLOW:
    //
    SEGEMEHL_ALIGN( reads, fasta, segemehl_index )

    segemehl_filter = SEGEMEHL_ALIGN.out.results.map{ meta, results ->  meta.tool = "segemehl"; return [ meta, results ] }

    SEGEMEHL_FILTER( segemehl_filter, bsj_reads )

    //
    // STAR WORFKLOW:
    //

    STAR_1ST_PASS( reads, star_index, gtf, true, '', '' )

    sjdb = STAR_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect()

    SJDB( sjdb, bsj_reads )

    STAR_2ND_PASS( reads, star_index, SJDB.out.sjtab, true, '', '' )

    //
    // CIRCEXPLORER2 WORKFLOW:
    //

    CIRCEXPLORER2_REFERENCE( gtf )

    CIRCEXPLORER2_PARSE( STAR_2ND_PASS.out.junction )

    CIRCEXPLORER2_ANNOTATE( CIRCEXPLORER2_PARSE.out.junction, fasta, CIRCEXPLORER2_REFERENCE.out.txt )

    circexplorer2_filter = CIRCEXPLORER2_ANNOTATE.out.txt.map{ meta, txt -> meta.tool = "circexplorer2"; return [ meta, txt ] }

    CIRCEXPLORER2_FILTER( circexplorer2_filter, bsj_reads )

    //
    // CIRCRNA_FINDER WORKFLOW:
    //

    CIRCRNA_FINDER_FILTER( STAR_2ND_PASS.out.bam.join( STAR_2ND_PASS.out.junction).join(STAR_2ND_PASS.out.tab), fasta, bsj_reads )

    //
    // ANNOTATION WORKFLOW:
    //

    circrna_filtered = CIRCEXPLORER2_FILTER.out.results.mix(SEGEMEHL_FILTER.out.results).view()

    ANNOTATION( circrna_filtered, gtf )

    // collect versions
    ch_versions = ch_versions.mix(CIRCEXPLORER2_REFERENCE.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_PARSE.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)
    ch_versions = ch_versions.mix(STAR_1ST_PASS.out.versions)

    emit:

    versions = ch_versions
}
