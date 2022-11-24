
include { ANNOTATION                       } from '../../modules/local/annotation/main'
include { BOWTIE2_ALIGN as FIND_CIRC_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_VIEW                    } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX                   } from '../../modules/nf-core/samtools/index/main'
include { FIND_CIRC_ANCHORS                } from '../../modules/local/find_circ/anchors/main'
include { FIND_CIRC                        } from '../../modules/local/find_circ/find_circ/main'
include { FIND_CIRC_FILTER                 } from '../../modules/local/find_circ/filter/main'
include { CIRIQUANT_YML                    } from '../../modules/local/ciriquant/yml/main'
include { CIRIQUANT                        } from '../../modules/local/ciriquant/ciriquant/main'
include { CIRIQUANT_FILTER                 } from '../../modules/local/ciriquant/filter/main'
include { CIRCEXPLORER2_REFERENCE as CIRCEXPLORER2_REFERENCE } from '../../modules/local/circexplorer2/reference/main'
include { CIRCEXPLORER2_PARSE as CIRCEXPLORER2_PARSE } from '../../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE as CIRCEXPLORER2_ANNOTATE } from '../../modules/nf-core/circexplorer2/annotate/main'
include { CIRCEXPLORER2_FILTER as CIRCEXPLORER2_FILTER } from '../../modules/local/circexplorer2/filter/main'
include { CIRCRNA_FINDER_FILTER            } from '../../modules/local/circrna_finder/filter/main'
include { SEGEMEHL_ALIGN                   } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_FILTER                  } from '../../modules/local/segemehl/filter/main'
include { STAR_ALIGN as STAR_1ST_PASS      } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_2ND_PASS      } from '../../modules/nf-core/star/align/main'
include { SJDB as SJDB                     } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE1_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE1_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE1_SJDB           } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE2_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE2_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE2_SJDB           } from '../../modules/local/star/sjdb/main'
include { DCC                              } from '../../modules/local/dcc/dcc/main'
include { DCC_FILTER                       } from '../../modules/local/dcc/filter/main'
include { CIRCEXPLORER2_REFERENCE as MAPSPLICE_REFERENCE } from '../../modules/local/circexplorer2/reference/main'
include { MAPSPLICE_ALIGN                  } from '../../modules/local/mapsplice/align/main'
include { CIRCEXPLORER2_PARSE as MAPSPLICE_PARSE } from '../../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE as MAPSPLICE_ANNOTATE } from '../../modules/nf-core/circexplorer2/annotate/main'
include { CIRCEXPLORER2_FILTER as MAPSPLICE_FILTER } from '../../modules/local/circexplorer2/filter/main'
include { FASTA                            } from '../../modules/local/fasta/main'
include { MERGE_TOOLS                      } from '../../modules/local/count_matrix/merge_tools/main'
include { COUNTS_COMBINED                  } from '../../modules/local/count_matrix/combined/main'
include { COUNTS_SINGLE                    } from '../../modules/local/count_matrix/single/main'

workflow CIRCRNA_DISCOVERY {

    take:
    reads
    fasta
    gtf
    bowtie_index
    bowtie2_index
    bwa_index
    chromosomes
    hisat2_index
    segemehl_index
    star_index
    bsj_reads
    tool_filter

    main:
    ch_versions = Channel.empty()

    //
    // SEGEMEHL WORKFLOW:
    //
    SEGEMEHL_ALIGN( reads, fasta, segemehl_index )
    segemehl_filter = SEGEMEHL_ALIGN.out.results.map{ meta, results ->  meta.tool = "segemehl"; return [ meta, results ] }
    SEGEMEHL_FILTER( segemehl_filter, bsj_reads )

    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)

    //
    // STAR WORFKLOW:
    //

    STAR_1ST_PASS( reads, star_index, gtf, true, '', '' )
    sjdb = STAR_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect()
    SJDB( sjdb, bsj_reads )
    STAR_2ND_PASS( reads, star_index, SJDB.out.sjtab, true, '', '' )

    ch_versions = ch_versions.mix(STAR_1ST_PASS.out.versions)

    //
    // CIRCEXPLORER2 WORKFLOW:
    //

    CIRCEXPLORER2_REFERENCE( gtf )
    CIRCEXPLORER2_PARSE( STAR_2ND_PASS.out.junction )
    CIRCEXPLORER2_ANNOTATE( CIRCEXPLORER2_PARSE.out.junction, fasta, CIRCEXPLORER2_REFERENCE.out.txt )
    circexplorer2_filter = CIRCEXPLORER2_ANNOTATE.out.txt.map{ meta, txt -> meta.tool = "circexplorer2"; return [ meta, txt ] }
    CIRCEXPLORER2_FILTER( circexplorer2_filter, bsj_reads )

    ch_versions = ch_versions.mix(CIRCEXPLORER2_REFERENCE.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_PARSE.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_ANNOTATE.out.versions)

    //
    // CIRCRNA_FINDER WORKFLOW:
    //

    circrna_finder_stage = STAR_2ND_PASS.out.sam.join( STAR_2ND_PASS.out.junction).join(STAR_2ND_PASS.out.tab)
    circrna_finder_filter = circrna_finder_stage.map{ meta, sam, junction, tab -> meta.tool = "circrna_finder"; return [ meta, sam, junction, tab ] }
    CIRCRNA_FINDER_FILTER( circrna_finder_filter, fasta, bsj_reads )

    ch_versions = ch_versions.mix(CIRCRNA_FINDER_FILTER.out.versions)

    //
    // FIND_CIRC WORKFLOW:
    //

    FIND_CIRC_ALIGN( reads, bowtie2_index.collect(), false, true )
    SAMTOOLS_INDEX( FIND_CIRC_ALIGN.out.bam )
    SAMTOOLS_VIEW( FIND_CIRC_ALIGN.out.bam.join( SAMTOOLS_INDEX.out.bai ), fasta, [] )
    FIND_CIRC_ANCHORS( SAMTOOLS_VIEW.out.bam )
    FIND_CIRC( FIND_CIRC_ANCHORS.out.anchors, bowtie2_index.collect(), fasta, chromosomes )
    find_circ_filter = FIND_CIRC.out.bed.map{ meta, bed -> meta.tool = "find_circ"; return [ meta, bed ] }
    FIND_CIRC_FILTER( find_circ_filter, bsj_reads )

    ch_versions = ch_versions.mix(FIND_CIRC_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_versions = ch_versions.mix(FIND_CIRC_ANCHORS.out.versions)
    ch_versions = ch_versions.mix(FIND_CIRC_FILTER.out.versions)

    //
    // CIRIQUANT WORKFLOW:
    //

    CIRIQUANT_YML( gtf, fasta, bwa_index.map{ meta, index -> return index }, hisat2_index )
    CIRIQUANT( reads, CIRIQUANT_YML.out.yml.collect() )
    CIRIQUANT_FILTER( CIRIQUANT.out.gtf.map{ meta, gtf -> meta.tool = "ciriquant"; return [ meta, gtf ] }, bsj_reads )

    ch_versions = ch_versions.mix(CIRIQUANT.out.versions)

    //
    // DCC WORKFLOW
    //

    mate1 = reads.map{ meta, reads -> return [ meta, reads[0] ] }
    DCC_MATE1_1ST_PASS( mate1, star_index, gtf, true, '', '' )
    DCC_MATE1_SJDB( DCC_MATE1_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect(), bsj_reads )
    DCC_MATE1_2ND_PASS( mate1, star_index, DCC_MATE1_SJDB.out.sjtab, true, '', '' )

    mate2 = reads.map{ meta, reads -> return [ meta, reads[1] ] }
    DCC_MATE2_1ST_PASS( mate2, star_index, gtf, true, '', '' )
    DCC_MATE2_SJDB( DCC_MATE2_1ST_PASS.out.tab.map{ meta, tab -> return [ tab ] }.collect(), bsj_reads )
    DCC_MATE2_2ND_PASS( mate2, star_index, DCC_MATE2_SJDB.out.sjtab, true, '', '' )

    dcc = STAR_2ND_PASS.out.junction.join( DCC_MATE1_2ND_PASS.out.junction ).join( DCC_MATE2_2ND_PASS.out.junction )
    DCC( dcc, fasta, gtf )
    DCC_FILTER( DCC.out.txt.map{ meta, txt -> meta.tool = "dcc"; return [ meta, txt ] }, bsj_reads )

    ch_versions = ch_versions.mix(DCC_MATE1_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC.out.versions)

    //
    // MAPSPLICE WORKFLOW:
    //

    MAPSPLICE_REFERENCE( gtf )
    MAPSPLICE_ALIGN( reads, bowtie_index.collect(), chromosomes, gtf )
    MAPSPLICE_PARSE( MAPSPLICE_ALIGN.out.raw_fusions )
    MAPSPLICE_ANNOTATE( MAPSPLICE_PARSE.out.junction, fasta, MAPSPLICE_REFERENCE.out.txt )
    mapsplice_filter = MAPSPLICE_ANNOTATE.out.txt.map{ meta, txt -> meta.tool = "mapsplice"; return [ meta, txt ] }
    MAPSPLICE_FILTER( mapsplice_filter, bsj_reads )

    ch_versions = ch_versions.mix(MAPSPLICE_REFERENCE.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_ALIGN.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_PARSE.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_ANNOTATE.out.versions)

    //
    // ANNOTATION WORKFLOW:
    //

    circrna_filtered = CIRCEXPLORER2_FILTER.out.results.mix(SEGEMEHL_FILTER.out.results, CIRCRNA_FINDER_FILTER.out.results, FIND_CIRC_FILTER.out.results, CIRIQUANT_FILTER.out.results, DCC_FILTER.out.results, MAPSPLICE_FILTER.out.results )
    ANNOTATION( circrna_filtered, gtf )

    ch_versions = ch_versions.mix(ANNOTATION.out.versions)

    //
    // FASTA WORKFLOW:
    //

    FASTA( ANNOTATION.out.bed, fasta )

    ch_versions = ch_versions.mix(FASTA.out.versions)

    //
    // COUNT MATRIX WORKFLOW:
    //

    ch_matrix = CIRCEXPLORER2_FILTER.out.matrix.mix(SEGEMEHL_FILTER.out.matrix, CIRCRNA_FINDER_FILTER.out.matrix, FIND_CIRC_FILTER.out.matrix, CIRIQUANT_FILTER.out.matrix, DCC_FILTER.out.matrix, MAPSPLICE_FILTER.out.matrix )
    tools_selected = params.tool.split(',').collect{it.trim().toLowerCase()}

    if( tools_selected.size() > 1){

        MERGE_TOOLS( ch_matrix.groupTuple(by:[0,0]), tool_filter )

        COUNTS_COMBINED( MERGE_TOOLS.out.merged.collect() )

    }else{

        ch_matrix.map{ meta, bed -> var = [:]; var.tool = meta.tool; return [ var, bed ] }.collect().groupTuple().view()
        COUNTS_SINGLE( ch_matrix.map{ meta, bed -> var = [:]; var.tool = meta.tool; return [ var, bed ] }.groupTuple() )

    }

    emit:
    circrna_bed12 = ANNOTATION.out.bed
    fasta = FASTA.out.analysis_fasta
    versions = ch_versions
}
