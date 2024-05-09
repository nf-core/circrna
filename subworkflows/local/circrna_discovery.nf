
include { ANNOTATION                                     } from '../../modules/local/annotation/full_annotation/main'
include { GNU_SORT as COMBINE_ANNOTATION_BEDS            } from '../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as COMBINE_ANNOTATION_GTFS            } from '../../modules/nf-core/gnu/sort/main'
include { GAWK as REMOVE_SCORE_STRAND                    } from '../../modules/nf-core/gawk/main'
include { BEDTOOLS_INTERSECT as INTERSECT_ANNOTATION     } from '../../modules/nf-core/bedtools/intersect/main'
include { BOWTIE2_ALIGN as FIND_CIRC_ALIGN               } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_VIEW                                  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX                                 } from '../../modules/nf-core/samtools/index/main'
include { FIND_CIRC_ANCHORS                              } from '../../modules/local/find_circ/anchors/main'
include { FIND_CIRC                                      } from '../../modules/local/find_circ/find_circ/main'
include { FIND_CIRC_FILTER                               } from '../../modules/local/find_circ/filter/main'
include { CIRIQUANT                                      } from '../../modules/local/ciriquant/ciriquant/main'
include { CIRIQUANT_FILTER                               } from '../../modules/local/ciriquant/filter/main'
include { CIRCRNA_FINDER_FILTER                          } from '../../modules/local/circrna_finder/filter/main'
include { SEGEMEHL_ALIGN                                 } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_FILTER                                } from '../../modules/local/segemehl/filter/main'
include { STAR_ALIGN as STAR_1ST_PASS                    } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_2ND_PASS                    } from '../../modules/nf-core/star/align/main'
include { SJDB as STAR_SJDB                              } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE1_1ST_PASS               } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE1_2ND_PASS               } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE1_SJDB                         } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE2_1ST_PASS               } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE2_2ND_PASS               } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE2_SJDB                         } from '../../modules/local/star/sjdb/main'
include { DCC                                            } from '../../modules/local/dcc/dcc/main'
include { DCC_FILTER                                     } from '../../modules/local/dcc/filter/main'
include { MAPSPLICE_ALIGN                                } from '../../modules/local/mapsplice/align/main'
include { FASTA                                          } from '../../modules/local/fasta/main'
include { MERGE_TOOLS                                    } from '../../modules/local/count_matrix/merge_tools/main'
include { COUNTS_COMBINED                                } from '../../modules/local/count_matrix/combined/main'
include { CIRCEXPLORER2_REFERENCE as CIRCEXPLORER2_REF   } from '../../modules/local/circexplorer2/reference/main'
include { CIRCEXPLORER2_PARSE as CIRCEXPLORER2_PAR       } from '../../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE as CIRCEXPLORER2_ANN    } from '../../modules/nf-core/circexplorer2/annotate/main'
include { CIRCEXPLORER2_FILTER as CIRCEXPLORER2_FLT      } from '../../modules/local/circexplorer2/filter/main'
include { CIRCEXPLORER2_REFERENCE as MAPSPLICE_REFERENCE } from '../../modules/local/circexplorer2/reference/main'
include { CIRCEXPLORER2_PARSE as MAPSPLICE_PARSE         } from '../../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE as MAPSPLICE_ANNOTATE   } from '../../modules/nf-core/circexplorer2/annotate/main'
include { CIRCEXPLORER2_FILTER as MAPSPLICE_FILTER       } from '../../modules/local/circexplorer2/filter/main'

workflow CIRCRNA_DISCOVERY {

    take:
    reads
    ch_fasta
    ch_gtf
    bowtie_index
    bowtie2_index
    bwa_index
    chromosomes
    hisat2_index
    segemehl_index
    star_index
    bsj_reads
    tool_filter
    duplicates_fun
    exon_boundary

    main:
    ch_versions = Channel.empty()
    fasta       = ch_fasta.map{meta, fasta -> fasta}
    gtf         = ch_gtf.map{meta, gtf -> gtf}

    //
    // SEGEMEHL WORKFLOW:
    //
    SEGEMEHL_ALIGN( reads, fasta, segemehl_index )
    segemehl_filter = SEGEMEHL_ALIGN.out.results.map{ meta, results ->  [ meta + [tool: "segemehl"], results ] }
    SEGEMEHL_FILTER( segemehl_filter, bsj_reads )

    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SEGEMEHL_FILTER.out.versions)

    //
    // STAR WORFKLOW:
    //

    // Define variables here, star_ignore_sjdbgtf not supposed to be toggled by user.
    star_ignore_sjdbgtf = true
    seq_center     = params.seq_center ?: ''
    seq_platform   = ''

    STAR_1ST_PASS( reads, star_index, ch_gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
    sjdb = STAR_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "star_sjdb"], it]}
    STAR_SJDB( sjdb, bsj_reads )
    STAR_2ND_PASS( reads, star_index, STAR_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )

    ch_versions = ch_versions.mix(STAR_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(STAR_SJDB.out.versions)
    ch_versions = ch_versions.mix(STAR_2ND_PASS.out.versions)

    //
    // CIRCEXPLORER2 WORKFLOW:
    //

    CIRCEXPLORER2_REF( gtf )
    CIRCEXPLORER2_PAR( STAR_2ND_PASS.out.junction )
    CIRCEXPLORER2_ANN( CIRCEXPLORER2_PAR.out.junction, fasta, CIRCEXPLORER2_REF.out.txt )
    circexplorer2_filter = CIRCEXPLORER2_ANN.out.txt.map{ meta, txt -> [ meta + [tool: "circexplorer2"], txt ] }
    CIRCEXPLORER2_FLT( circexplorer2_filter, bsj_reads )

    ch_versions = ch_versions.mix(CIRCEXPLORER2_REF.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_PAR.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_ANN.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_FLT.out.versions)

    //
    // CIRCRNA_FINDER WORKFLOW:
    //

    circrna_finder_stage = STAR_2ND_PASS.out.sam.join( STAR_2ND_PASS.out.junction).join(STAR_2ND_PASS.out.tab)
    circrna_finder_filter = circrna_finder_stage.map{ meta, sam, junction, tab -> [ meta + [tool: "circrna_finder"], sam, junction, tab ] }
    CIRCRNA_FINDER_FILTER( circrna_finder_filter, fasta, bsj_reads )

    ch_versions = ch_versions.mix(CIRCRNA_FINDER_FILTER.out.versions)

    //
    // FIND_CIRC WORKFLOW:
    //

    FIND_CIRC_ALIGN( reads, bowtie2_index, ch_fasta, false, true )
    SAMTOOLS_INDEX( FIND_CIRC_ALIGN.out.bam )
    SAMTOOLS_VIEW( FIND_CIRC_ALIGN.out.bam.join( SAMTOOLS_INDEX.out.bai ), ch_fasta, [] )
    FIND_CIRC_ANCHORS( SAMTOOLS_VIEW.out.bam )
    FIND_CIRC( FIND_CIRC_ANCHORS.out.anchors, bowtie2_index, fasta )
    find_circ_filter = FIND_CIRC.out.bed.map{ meta, bed -> [ meta + [tool: "find_circ"], bed ] }
    FIND_CIRC_FILTER( find_circ_filter, bsj_reads )

    ch_versions = ch_versions.mix(FIND_CIRC_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_versions = ch_versions.mix(FIND_CIRC_ANCHORS.out.versions)
    ch_versions = ch_versions.mix(FIND_CIRC.out.versions)
    ch_versions = ch_versions.mix(FIND_CIRC_FILTER.out.versions)

    //
    // CIRIQUANT WORKFLOW:
    //

    // only need path to bwa, only need path to hisat2.
    // do not want to upset the collect declr for all indices just for this.
    CIRIQUANT( reads, ch_gtf, ch_fasta, bwa_index, hisat2_index )
    CIRIQUANT_FILTER( CIRIQUANT.out.gtf.map{ meta, gtf -> [ meta + [tool: "ciriquant"], gtf ] }, bsj_reads )

    ch_versions = ch_versions.mix(CIRIQUANT.out.versions)
    ch_versions = ch_versions.mix(CIRIQUANT_FILTER.out.versions)

    //
    // DCC WORKFLOW
    //

    mate1 = reads.filter{ meta, reads -> !meta.single_end }.map{ meta, reads -> return [ [id: meta.id, single_end: true], reads[0] ] }
    DCC_MATE1_1ST_PASS( mate1, star_index, ch_gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    DCC_MATE1_SJDB( DCC_MATE1_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "mate1_sjdb"], it]}, bsj_reads )
    DCC_MATE1_2ND_PASS( mate1, star_index, DCC_MATE1_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )

    mate2 = reads.filter{ meta, reads -> !meta.single_end }.map{ meta, reads -> return [ [id: meta.id, single_end: true], reads[1] ] }
    DCC_MATE2_1ST_PASS( mate2, star_index, ch_gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    DCC_MATE2_SJDB( DCC_MATE2_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "mate2_sjdb"], it]}, bsj_reads )
    DCC_MATE2_2ND_PASS( mate2, star_index, DCC_MATE2_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )

    dcc_stage = STAR_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, meta, junction]}
        .join(
            DCC_MATE1_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, junction] },
            remainder: true
        )
        .join(
            DCC_MATE2_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, junction] },
            remainder: true
        )
        .map{ id, meta, junction, mate1, mate2 -> return [ meta, junction, mate1, mate2 ]}

    dcc = dcc_stage.map{ it ->  [ it[0], it[1], it[2] ?: [], it[3] ?: [] ] }
    DCC( dcc, fasta, gtf )
    DCC_FILTER( DCC.out.txt.map{ meta, txt -> [ meta + [tool: "dcc"], txt ] }, bsj_reads )

    ch_versions = ch_versions.mix(DCC_MATE1_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE1_SJDB.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE1_2ND_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE2_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE2_SJDB.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE2_2ND_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC.out.versions)
    ch_versions = ch_versions.mix(DCC_FILTER.out.versions)

    //
    // MAPSPLICE WORKFLOW:
    //

    MAPSPLICE_REFERENCE( gtf )
    MAPSPLICE_ALIGN( reads, bowtie_index, chromosomes, gtf )
    MAPSPLICE_PARSE( MAPSPLICE_ALIGN.out.raw_fusions )
    MAPSPLICE_ANNOTATE( MAPSPLICE_PARSE.out.junction, fasta, MAPSPLICE_REFERENCE.out.txt )
    mapsplice_filter = MAPSPLICE_ANNOTATE.out.txt.map{ meta, txt -> [ meta + [tool: "mapsplice"], txt ] }
    MAPSPLICE_FILTER( mapsplice_filter, bsj_reads )

    ch_versions = ch_versions.mix(MAPSPLICE_REFERENCE.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_ALIGN.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_PARSE.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(MAPSPLICE_FILTER.out.versions)

    //
    // ANNOTATION WORKFLOW:
    //

    ch_biotypes = Channel.fromPath("${projectDir}/bin/unwanted_biotypes.txt")

    circrna_filtered = CIRCEXPLORER2_FLT.out.results.mix(SEGEMEHL_FILTER.out.results,
                                                            CIRCRNA_FINDER_FILTER.out.results,
                                                            FIND_CIRC_FILTER.out.results,
                                                            CIRIQUANT_FILTER.out.results,
                                                            DCC_FILTER.out.results,
                                                            MAPSPLICE_FILTER.out.results)

    INTERSECT_ANNOTATION( circrna_filtered.combine(gtf), [[], []])
    ANNOTATION( INTERSECT_ANNOTATION.out.intersect, exon_boundary )
    COMBINE_ANNOTATION_BEDS(ANNOTATION.out.bed.map{ meta, bed -> bed}.collect().map{[[id: "annotation"], it]})
    REMOVE_SCORE_STRAND( COMBINE_ANNOTATION_BEDS.out.sorted, [])
    COMBINE_ANNOTATION_GTFS(ANNOTATION.out.gtf.map{ meta, gtf -> gtf}.collect().map{[[id: "annotation"], it]})

    ch_versions = ch_versions.mix(INTERSECT_ANNOTATION.out.versions)
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)
    ch_versions = ch_versions.mix(COMBINE_ANNOTATION_BEDS.out.versions)
    ch_versions = ch_versions.mix(REMOVE_SCORE_STRAND.out.versions)
    ch_versions = ch_versions.mix(COMBINE_ANNOTATION_GTFS.out.versions)

    //
    // FASTA WORKFLOW:
    //

    FASTA( ANNOTATION.out.bed, fasta )

    ch_versions = ch_versions.mix(FASTA.out.versions)

    //
    // COUNT MATRIX WORKFLOW:
    //

    ch_matrix = CIRCEXPLORER2_FLT.out.matrix.mix(SEGEMEHL_FILTER.out.matrix,
                                                    CIRCRNA_FINDER_FILTER.out.matrix,
                                                    FIND_CIRC_FILTER.out.matrix,
                                                    CIRIQUANT_FILTER.out.matrix,
                                                    DCC_FILTER.out.matrix,
                                                    MAPSPLICE_FILTER.out.matrix)

    tools_selected = params.tool.split(',').collect{it.trim().toLowerCase()}

    MERGE_TOOLS( ch_matrix.map{ meta, bed -> [ [id: meta.id], bed ] }.groupTuple(),
                tools_selected.size() > 1 ? tool_filter : 1, duplicates_fun )

    COUNTS_COMBINED( MERGE_TOOLS.out.merged.map{ meta, bed -> bed }.collect() )

    counts_bed = COUNTS_COMBINED.out.counts_bed
    counts_tsv = COUNTS_COMBINED.out.counts_tsv
    ch_versions = ch_versions.mix(MERGE_TOOLS.out.versions)
    ch_versions = ch_versions.mix(COUNTS_COMBINED.out.versions)

    emit:
    circrna_bed12 = ANNOTATION.out.bed
    fasta = FASTA.out.analysis_fasta
    annotation_bed = REMOVE_SCORE_STRAND.out.output
    annotation_gtf = COMBINE_ANNOTATION_GTFS.out.sorted
    counts_bed
    counts_tsv

    versions = ch_versions
}
