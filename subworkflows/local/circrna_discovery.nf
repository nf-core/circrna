// MODULES
include { ANNOTATION                                     } from '../../modules/local/annotation/full_annotation'
include { GNU_SORT as COMBINE_ANNOTATION_BEDS            } from '../../modules/nf-core/gnu/sort'
include { GNU_SORT as COMBINE_ANNOTATION_GTFS            } from '../../modules/nf-core/gnu/sort'
include { GAWK as REMOVE_SCORE_STRAND                    } from '../../modules/nf-core/gawk'
include { BEDTOOLS_INTERSECT as INTERSECT_ANNOTATION     } from '../../modules/nf-core/bedtools/intersect'
include { BOWTIE2_ALIGN as FIND_CIRC_ALIGN               } from '../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_VIEW                                  } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX                                 } from '../../modules/nf-core/samtools/index'
include { FIND_CIRC_ANCHORS                              } from '../../modules/local/find_circ/anchors'
include { FIND_CIRC                                      } from '../../modules/local/find_circ/find_circ'
include { FIND_CIRC_FILTER                               } from '../../modules/local/find_circ/filter'
include { CIRIQUANT                                      } from '../../modules/local/ciriquant/ciriquant'
include { CIRIQUANT_FILTER                               } from '../../modules/local/ciriquant/filter'
include { STAR_ALIGN as DCC_MATE1_1ST_PASS               } from '../../modules/nf-core/star/align'
include { STAR_ALIGN as DCC_MATE1_2ND_PASS               } from '../../modules/nf-core/star/align'
include { SJDB as DCC_MATE1_SJDB                         } from '../../modules/local/star/sjdb'
include { STAR_ALIGN as DCC_MATE2_1ST_PASS               } from '../../modules/nf-core/star/align'
include { STAR_ALIGN as DCC_MATE2_2ND_PASS               } from '../../modules/nf-core/star/align'
include { SJDB as DCC_MATE2_SJDB                         } from '../../modules/local/star/sjdb'
include { DCC                                            } from '../../modules/local/dcc/dcc'
include { DCC_FILTER                                     } from '../../modules/local/dcc/filter'
include { MAPSPLICE_ALIGN                                } from '../../modules/local/mapsplice/align'
include { MERGE_TOOLS                                    } from '../../modules/local/count_matrix/merge_tools'
include { COUNTS_COMBINED                                } from '../../modules/local/count_matrix/combined'
include { CIRCEXPLORER2_REFERENCE as MAPSPLICE_REFERENCE } from '../../modules/local/circexplorer2/reference'
include { CIRCEXPLORER2_PARSE as MAPSPLICE_PARSE         } from '../../modules/nf-core/circexplorer2/parse'
include { CIRCEXPLORER2_ANNOTATE as MAPSPLICE_ANNOTATE   } from '../../modules/nf-core/circexplorer2/annotate'
include { CIRCEXPLORER2_FILTER as MAPSPLICE_FILTER       } from '../../modules/local/circexplorer2/filter'
include { UPSET as UPSET_SAMPLES                         } from '../../modules/local/upset'
include { UPSET as UPSET_ALL                             } from '../../modules/local/upset'
include { BEDTOOLS_GETFASTA                              } from '../../modules/nf-core/bedtools/getfasta'
include { GAWK as ADD_BACKSPLICE                         } from '../../modules/nf-core/gawk'

// SUBWORKFLOWS
include { SEGEMEHL       } from './discovery/segemehl'
include { STAR2PASS      } from './discovery/star2pass'
include { CIRCEXPLORER2  } from './discovery/circexplorer2'
include { CIRCRNA_FINDER } from './discovery/circrna_finder'

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
    star_index
    bsj_reads
    tool_filter
    duplicates_fun
    exon_boundary

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    fasta       = ch_fasta.map{meta, fasta -> fasta}
    gtf         = ch_gtf.map{meta, gtf -> gtf}

    // STAR CONSTANT PARAMETERS
    star_ignore_sjdbgtf = true
    seq_center = params.seq_center ?: ''
    seq_platform = ''

    //
    // DISCOVERY TOOLS:
    //
    tools = params.tool.split(',').collect{it.trim().toLowerCase()}
    SEGEMEHL( reads, fasta, params.segemehl, bsj_reads )
    STAR2PASS( reads, star_index, ch_gtf, bsj_reads, star_ignore_sjdbgtf, seq_center, seq_platform )
    CIRCEXPLORER2( gtf, fasta, STAR2PASS.out.junction, bsj_reads )
    CIRCRNA_FINDER( fasta, STAR2PASS.out.sam, STAR2PASS.out.junction, STAR2PASS.out.tab, bsj_reads )

    ch_versions = ch_versions.mix(SEGEMEHL.out.versions)
    ch_versions = ch_versions.mix(STAR2PASS.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2.out.versions)
    ch_versions = ch_versions.mix(CIRCRNA_FINDER.out.versions)


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

    dcc_stage = STAR2PASS.out.junction.map{ meta, junction -> return [ meta.id, meta, junction]}
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
    // COUNT MATRIX WORKFLOW:
    //

    ch_matrix = CIRCEXPLORER2.out.matrix.mix(SEGEMEHL.out.matrix,
                                                    CIRCRNA_FINDER.out.matrix,
                                                    FIND_CIRC_FILTER.out.matrix,
                                                    CIRIQUANT_FILTER.out.matrix,
                                                    DCC_FILTER.out.matrix,
                                                    MAPSPLICE_FILTER.out.matrix)

    tools_selected = params.tool.split(',').collect{it.trim().toLowerCase()}

    MERGE_TOOLS( ch_matrix.map{ meta, bed -> [ [id: meta.id], bed ] }.groupTuple(),
                tools_selected.size() > 1 ? tool_filter : 1, duplicates_fun )
    COUNTS_COMBINED( MERGE_TOOLS.out.merged.map{ meta, bed -> bed }.collect() )

    ch_versions = ch_versions.mix(MERGE_TOOLS.out.versions)
    ch_versions = ch_versions.mix(COUNTS_COMBINED.out.versions)

    //
    // ANNOTATION WORKFLOW:
    //

    ch_biotypes = Channel.fromPath("${projectDir}/bin/unwanted_biotypes.txt")

    circrna_tools = CIRCEXPLORER2.out.results.mix(SEGEMEHL.out.results,
                                                            CIRCRNA_FINDER.out.results,
                                                            FIND_CIRC_FILTER.out.results,
                                                            CIRIQUANT_FILTER.out.results,
                                                            DCC_FILTER.out.results,
                                                            MAPSPLICE_FILTER.out.results)

    UPSET_SAMPLES( circrna_tools.map{ meta, bed -> [meta.id, meta.tool, bed]}
        .groupTuple()
        .map{ sample, tools, beds -> [[id: sample], tools, beds]} )
    UPSET_ALL( circrna_tools.map{ meta, bed -> ["all", meta.tool, bed] }
        .groupTuple()
        .map{ sample, tools, beds -> [[id: sample], tools, beds]} )

    ch_multiqc_files = ch_multiqc_files.mix(UPSET_SAMPLES.out.multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(UPSET_ALL.out.multiqc)
    ch_versions = ch_versions.mix(UPSET_SAMPLES.out.versions)
    ch_versions = ch_versions.mix(UPSET_ALL.out.versions)

    circrna_incl_merged = circrna_tools.mix(
        MERGE_TOOLS.out.merged.map{ meta, bed -> [meta + [tool: "merged"], bed] })

    INTERSECT_ANNOTATION( circrna_incl_merged.combine(gtf), [[], []])
    ANNOTATION( INTERSECT_ANNOTATION.out.intersect, exon_boundary )

    ch_annotation_bed_merged = ANNOTATION.out.bed.filter{ meta, bed -> meta.tool == "merged" }
    ch_annotation_gtf_merged = ANNOTATION.out.gtf.filter{ meta, gtf -> meta.tool == "merged" }

    COMBINE_ANNOTATION_BEDS(ch_annotation_bed_merged.map{ meta, bed -> bed}.collect().map{[[id: "annotation"], it]})
    REMOVE_SCORE_STRAND( COMBINE_ANNOTATION_BEDS.out.sorted, [])
    COMBINE_ANNOTATION_GTFS(ch_annotation_gtf_merged.map{ meta, gtf -> gtf}.collect().map{[[id: "annotation"], it]})

    ch_versions = ch_versions.mix(INTERSECT_ANNOTATION.out.versions)
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)
    ch_versions = ch_versions.mix(COMBINE_ANNOTATION_BEDS.out.versions)
    ch_versions = ch_versions.mix(REMOVE_SCORE_STRAND.out.versions)
    ch_versions = ch_versions.mix(COMBINE_ANNOTATION_GTFS.out.versions)

    //
    // FASTA WORKFLOW:
    //

    BEDTOOLS_GETFASTA( ch_annotation_bed_merged, fasta )
    ADD_BACKSPLICE( BEDTOOLS_GETFASTA.out.fasta, [])

    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)
    ch_versions = ch_versions.mix(ADD_BACKSPLICE.out.versions)

    emit:
    circrna_bed12 = ch_annotation_bed_merged
    fasta = ADD_BACKSPLICE.out.output
    annotation_bed = REMOVE_SCORE_STRAND.out.output
    annotation_gtf = COMBINE_ANNOTATION_GTFS.out.sorted
    counts_bed = COUNTS_COMBINED.out.counts_bed
    counts_tsv = COUNTS_COMBINED.out.counts_tsv

    multiqc_files = ch_multiqc_files
    versions = ch_versions
}
