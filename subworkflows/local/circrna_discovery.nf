// MODULES
include { GAWK as FILTER_BSJS    } from '../../modules/nf-core/gawk'
include { MERGE_TOOLS            } from '../../modules/local/count_matrix/merge_tools'
include { MERGE_SAMPLES          } from '../../modules/local/count_matrix/merge_samples'
include { UPSET as UPSET_SAMPLES } from '../../modules/local/upset'
include { UPSET as UPSET_ALL     } from '../../modules/local/upset'
include { BEDTOOLS_GETFASTA      } from '../../modules/nf-core/bedtools/getfasta'
include { GAWK as ADD_BACKSPLICE } from '../../modules/nf-core/gawk'

// SUBWORKFLOWS
include { SEGEMEHL       } from './discovery/segemehl'
include { STAR2PASS      } from './discovery/star2pass'
include { CIRCEXPLORER2  } from './discovery/circexplorer2'
include { CIRCRNA_FINDER } from './discovery/circrna_finder'
include { FIND_CIRC      } from './discovery/find_circ'
include { CIRIQUANT      } from './discovery/ciriquant'
include { DCC            } from './discovery/dcc'
include { MAPSPLICE      } from './discovery/mapsplice'
include { ANNOTATION     } from './discovery/annotation'

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
    ch_annotation
    bsj_reads
    tool_filter
    duplicates_fun
    exon_boundary

    main:
    ch_versions      = Channel.empty()
    ch_bed           = Channel.empty()
    ch_multiqc_files = Channel.empty()
    fasta            = ch_fasta.map{meta, fasta -> fasta}
    gtf              = ch_gtf.map{meta, gtf -> gtf}

    // STAR 2-PASS-MODE
    star_ignore_sjdbgtf = true
    seq_center = params.seq_center ?: ''
    seq_platform = ''
    STAR2PASS( reads, star_index, ch_gtf, bsj_reads, star_ignore_sjdbgtf, seq_center, seq_platform )
    ch_versions = ch_versions.mix(STAR2PASS.out.versions)

    //
    // DISCOVERY TOOLS:
    //
    tools_selected = params.tool.split(',').collect{it.trim().toLowerCase()}

    if (tools_selected.size() == 0) {
        error 'No tools selected for circRNA discovery.'
    }

    if (tools_selected.contains('segemehl')) {
        SEGEMEHL( reads, fasta, params.segemehl )
        ch_versions = ch_versions.mix(SEGEMEHL.out.versions)
        ch_bed      = ch_bed     .mix(SEGEMEHL.out.bed)
    }

    if (tools_selected.contains('circexplorer2')) {
        CIRCEXPLORER2( gtf, fasta, STAR2PASS.out.junction )
        ch_versions = ch_versions.mix(CIRCEXPLORER2.out.versions)
        ch_bed      = ch_bed     .mix(CIRCEXPLORER2.out.bed)
    }

    if (tools_selected.contains('circrna_finder')) {
        CIRCRNA_FINDER( fasta, STAR2PASS.out.sam, STAR2PASS.out.junction,
            STAR2PASS.out.tab )
        ch_versions = ch_versions.mix(CIRCRNA_FINDER.out.versions)
        ch_bed      = ch_bed     .mix(CIRCRNA_FINDER.out.bed)
    }

    if (tools_selected.contains('find_circ')) {
        FIND_CIRC( reads, bowtie2_index, ch_fasta )
        ch_versions = ch_versions.mix(FIND_CIRC.out.versions)
        ch_bed      = ch_bed     .mix(FIND_CIRC.out.bed)
    }

    if (tools_selected.contains('ciriquant')) {
        CIRIQUANT( reads, ch_gtf, ch_fasta, bwa_index, hisat2_index )
        ch_versions = ch_versions.mix(CIRIQUANT.out.versions)
        ch_bed      = ch_bed     .mix(CIRIQUANT.out.bed)
    }

    if (tools_selected.contains('dcc')) {
        DCC( reads, ch_fasta, ch_gtf, star_index, STAR2PASS.out.junction,
            star_ignore_sjdbgtf, seq_platform, seq_center, bsj_reads )
        ch_versions = ch_versions.mix(DCC.out.versions)
        ch_bed      = ch_bed     .mix(DCC.out.bed)
    }

    if (tools_selected.contains('mapsplice')) {
        MAPSPLICE( reads, gtf, fasta, bowtie_index, chromosomes,
            STAR2PASS.out.junction )
        ch_versions = ch_versions.mix(MAPSPLICE.out.versions)
        ch_bed      = ch_bed     .mix(MAPSPLICE.out.bed)
    }

    ch_bed = FILTER_BSJS( ch_bed, [] ).output
    ch_versions = ch_versions.mix(FILTER_BSJS.out.versions)

    //
    // CREATE COUNT MATRIX
    //

    MERGE_TOOLS( ch_bed.map{ meta, bed -> [ [id: meta.id], bed ] }.groupTuple(),
                tools_selected.size() > 1 ? tool_filter : 1, duplicates_fun )
    MERGE_SAMPLES( MERGE_TOOLS.out.merged.map{ meta, bed -> bed }.collect() )

    ch_bed_incl_merged = ch_bed.mix(
        MERGE_TOOLS.out.merged.map{ meta, bed -> [meta + [tool: "merged"], bed] })

    ch_versions = ch_versions.mix(MERGE_TOOLS.out.versions)
    ch_versions = ch_versions.mix(MERGE_SAMPLES.out.versions)

    //
    // UPSET PLOTS
    //

    UPSET_SAMPLES( ch_bed.map{ meta, bed -> [meta.id, meta.tool, bed]}
        .groupTuple()
        .map{ sample, tools, beds -> [[id: sample], tools, beds]} )
    UPSET_ALL( ch_bed.map{ meta, bed -> ["all", meta.tool, bed] }
        .groupTuple()
        .map{ sample, tools, beds -> [[id: sample], tools, beds]} )

    ch_multiqc_files = ch_multiqc_files.mix(UPSET_SAMPLES.out.multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(UPSET_ALL.out.multiqc)
    ch_versions = ch_versions.mix(UPSET_SAMPLES.out.versions)
    ch_versions = ch_versions.mix(UPSET_ALL.out.versions)

    //
    // ANNOTATION WORKFLOW:
    //

    ANNOTATION( ch_bed_incl_merged, gtf, exon_boundary, ch_annotation )
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)

    //
    // FASTA WORKFLOW:
    //

    BEDTOOLS_GETFASTA( ANNOTATION.out.merged_bed, fasta )
    ADD_BACKSPLICE( BEDTOOLS_GETFASTA.out.fasta, [])

    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)
    ch_versions = ch_versions.mix(ADD_BACKSPLICE.out.versions)

    emit:
    circrna_bed12  = ANNOTATION.out.merged_bed
    fasta          = ADD_BACKSPLICE.out.output
    annotation_bed = ANNOTATION.out.bed
    annotation_gtf = ANNOTATION.out.gtf
    counts_bed     = MERGE_SAMPLES.out.counts_bed
    counts_tsv     = MERGE_SAMPLES.out.counts_tsv

    multiqc_files  = ch_multiqc_files
    versions       = ch_versions
}
