// MODULES
include { GAWK as FILTER_BSJS                 } from '../../modules/nf-core/gawk'
include { GAWK as MASK_SCORES                 } from '../../modules/nf-core/gawk'
include { GNU_SORT as CONCAT_TOOLS_PER_SAMPLE } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_GROUPBY as COUNT_TOOLS     } from '../../modules/nf-core/bedtools/groupby'
include { GAWK as FILTER_MIN_TOOLS            } from '../../modules/nf-core/gawk'
include { GNU_SORT as CONCAT_SAMPLES          } from '../../modules/nf-core/gnu/sort'
include { UPSET as UPSET_SAMPLES              } from '../../modules/local/upset'
include { UPSET as UPSET_ALL                  } from '../../modules/local/upset'
include { BEDTOOLS_GETFASTA                   } from '../../modules/nf-core/bedtools/getfasta'

// SUBWORKFLOWS
include { SEGEMEHL       } from './discovery/segemehl'
include { STAR2PASS      } from './discovery/star2pass'
include { CIRCEXPLORER2  } from './discovery/circexplorer2'
include { CIRCRNA_FINDER } from './discovery/circrna_finder'
include { FIND_CIRC      } from './discovery/find_circ'
include { CIRIQUANT      } from './discovery/ciriquant'
include { DCC            } from './discovery/dcc'
include { MAPSPLICE      } from './discovery/mapsplice'

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
    // MERGE BED FILES
    //

    MASK_SCORES( ch_bed, [] )
    ch_versions = ch_versions.mix(MASK_SCORES.out.versions)
    ch_bsj_bed_per_sample_tool = MASK_SCORES.out.output

    CONCAT_TOOLS_PER_SAMPLE(
        MASK_SCORES.out.output.map{ meta, bed -> [ [id: meta.id], bed ] }.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCAT_TOOLS_PER_SAMPLE.out.versions)

    COUNT_TOOLS( CONCAT_TOOLS_PER_SAMPLE.out.sorted, 5 )
    ch_versions = ch_versions.mix(COUNT_TOOLS.out.versions)

    FILTER_MIN_TOOLS( COUNT_TOOLS.out.bed, [] )
    ch_versions = ch_versions.mix(FILTER_MIN_TOOLS.out.versions)
    ch_bsj_bed_per_sample = FILTER_MIN_TOOLS.out.output

    CONCAT_SAMPLES(
        ch_bsj_bed_per_sample.map{ meta, bed -> [[id: "all"], bed] }.groupTuple()
    )
    ch_versions = ch_versions.mix(CONCAT_SAMPLES.out.versions)
    ch_bsj_bed_combined = CONCAT_SAMPLES.out.sorted

    //
    // UPSET PLOTS
    //

    UPSET_SAMPLES( ch_bed.map{ meta, bed -> [meta.id, meta.tool, bed]}
        .groupTuple()
        .map{ sample, tools, beds -> [[id: sample], tools, beds]} )
    ch_multiqc_files = ch_multiqc_files.mix(UPSET_SAMPLES.out.multiqc)
    ch_versions = ch_versions.mix(UPSET_SAMPLES.out.versions)

    UPSET_ALL( ch_bed.map{ meta, bed -> ["all", meta.tool, bed] }
        .groupTuple()
        .map{ sample, tools, beds -> [[id: sample], tools, beds]} )
    ch_multiqc_files = ch_multiqc_files.mix(UPSET_ALL.out.multiqc)
    ch_versions = ch_versions.mix(UPSET_ALL.out.versions)

    //
    // FASTA WORKFLOW:
    //

    BEDTOOLS_GETFASTA(
        ch_bsj_bed_combined.mix(ch_bsj_bed_per_sample).mix(ch_bsj_bed_per_sample_tool),
        fasta
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

    ch_bsj_fasta = BEDTOOLS_GETFASTA.out.fasta.branch{
        meta, fasta ->
            combined: meta.id == 'all'
            per_sample_tool: meta.containsKey('tool')
            per_sample: true
    }

    emit:
    bsj_bed_combined          = ch_bsj_bed_combined
    bsj_fasta_combined        = ch_bsj_fasta.combined

    bsj_bed_per_sample        = ch_bsj_bed_per_sample
    bsj_fasta_per_sample      = ch_bsj_fasta.per_sample

    bsj_bed_per_sample_tool   = ch_bsj_bed_per_sample_tool
    bsj_fasta_per_sample_tool = ch_bsj_fasta.per_sample_tool

    multiqc_files             = ch_multiqc_files
    versions                  = ch_versions
}
