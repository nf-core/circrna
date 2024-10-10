// MODULES
include { GAWK as EXTRACT_COUNTS                         } from '../../modules/nf-core/gawk'
include { CSVTK_JOIN as COMBINE_COUNTS_PER_TOOL          } from '../../modules/nf-core/csvtk/join'
include { GAWK as FILTER_BSJS                            } from '../../modules/nf-core/gawk'
include { GAWK as BED_ADD_SAMPLE_TOOL                    } from '../../modules/nf-core/gawk'
include { COMBINEBEDS_FILTER as COMBINE_TOOLS_PER_SAMPLE } from '../../modules/local/combinebeds/filter'
include { COMBINEBEDS_FILTER as COMBINE_SAMPLES          } from '../../modules/local/combinebeds/filter'
include { AGAT_SPADDINTRONS as ADD_INTRONS          } from '../../modules/nf-core/agat/spaddintrons'
include { GAWK as EXTRACT_EXONS_INTRONS             } from '../../modules/nf-core/gawk'
include { BEDTOOLS_GETFASTA as FASTA_COMBINED            } from '../../modules/nf-core/bedtools/getfasta'
include { BEDTOOLS_GETFASTA as FASTA_PER_SAMPLE          } from '../../modules/nf-core/bedtools/getfasta'
include { BEDTOOLS_GETFASTA as FASTA_PER_SAMPLE_TOOL     } from '../../modules/nf-core/bedtools/getfasta'
include { FAIL_ON_EMPTY                                  } from '../../modules/local/fail_on_empty'

// SUBWORKFLOWS
include { SEGEMEHL                               } from './detection_tools/segemehl'
include { STAR2PASS                              } from './detection_tools/star2pass'
include { CIRCEXPLORER2                          } from './detection_tools/circexplorer2'
include { CIRCRNA_FINDER                         } from './detection_tools/circrna_finder'
include { FIND_CIRC                              } from './detection_tools/find_circ'
include { CIRIQUANT                              } from './detection_tools/ciriquant'
include { DCC                                    } from './detection_tools/dcc'
include { MAPSPLICE                              } from './detection_tools/mapsplice'
include { ANNOTATION as ANNOTATE_COMBINED        } from './annotation'
include { ANNOTATION as ANNOTATE_PER_SAMPLE      } from './annotation'
include { ANNOTATION as ANNOTATE_PER_SAMPLE_TOOL } from './annotation'

workflow BSJ_DETECTION {

    take:
    reads
    ch_fasta
    ch_gtf
    ch_annotation
    bowtie_index
    bowtie2_index
    bwa_index
    chromosomes
    hisat2_index
    star_index
    bsj_reads

    main:
    ch_versions                = Channel.empty()
    ch_bsj_bed_per_sample_tool = Channel.empty()
    ch_multiqc_files           = Channel.empty()
    fasta                      = ch_fasta.map{meta, fasta -> fasta}
    gtf                        = ch_gtf.map{meta, gtf -> gtf}

    // STAR 2-PASS-MODE
    star_ignore_sjdbgtf = true
    seq_center = params.seq_center ?: ''
    seq_platform = ''
    STAR2PASS( reads, star_index, ch_gtf, bsj_reads, star_ignore_sjdbgtf, seq_center, seq_platform )
    ch_versions = ch_versions.mix(STAR2PASS.out.versions)

    //
    // DISCOVERY TOOLS:
    //
    tools_selected = params.tools.split(',').collect{it.trim().toLowerCase()}

    if (tools_selected.size() == 0) {
        error 'No tools selected for circRNA discovery.'
    }

    if (tools_selected.contains('segemehl')) {
        SEGEMEHL( reads, fasta, params.segemehl )
        ch_versions                = ch_versions.mix(SEGEMEHL.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(SEGEMEHL.out.bed)
    }

    if (tools_selected.contains('circexplorer2')) {
        CIRCEXPLORER2( gtf, fasta, STAR2PASS.out.junction )
        ch_versions                = ch_versions.mix(CIRCEXPLORER2.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCEXPLORER2.out.bed)
    }

    if (tools_selected.contains('circrna_finder')) {
        CIRCRNA_FINDER( fasta, STAR2PASS.out.sam, STAR2PASS.out.junction,
            STAR2PASS.out.tab )
        ch_versions                = ch_versions.mix(CIRCRNA_FINDER.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCRNA_FINDER.out.bed)
    }

    if (tools_selected.contains('find_circ')) {
        FIND_CIRC( reads, bowtie2_index, ch_fasta )
        ch_versions                = ch_versions.mix(FIND_CIRC.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(FIND_CIRC.out.bed)
    }

    if (tools_selected.contains('ciriquant')) {
        CIRIQUANT( reads, ch_gtf, ch_fasta, bwa_index, hisat2_index )
        ch_versions                = ch_versions.mix(CIRIQUANT.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRIQUANT.out.bed)
    }

    if (tools_selected.contains('dcc')) {
        DCC( reads, ch_fasta, ch_gtf, star_index, STAR2PASS.out.junction,
            star_ignore_sjdbgtf, seq_platform, seq_center, bsj_reads )
        ch_versions                = ch_versions.mix(DCC.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(DCC.out.bed)
    }

    if (tools_selected.contains('mapsplice')) {
        MAPSPLICE( reads, gtf, fasta, bowtie_index, chromosomes,
            STAR2PASS.out.junction )
        ch_versions                = ch_versions.mix(MAPSPLICE.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(MAPSPLICE.out.bed)
    }

    //
    // QUANTIFY BSJs PER TOOL
    //

    EXTRACT_COUNTS( ch_bsj_bed_per_sample_tool, [] )
    ch_versions = ch_versions.mix(EXTRACT_COUNTS.out.versions)

    COMBINE_COUNTS_PER_TOOL( EXTRACT_COUNTS.out.output
        .map{ meta, bed -> [[id: meta.tool], bed]}
        .groupTuple() )
    ch_versions = ch_versions.mix(COMBINE_COUNTS_PER_TOOL.out.versions)

    //
    // APPLY bsj_reads FILTER
    //

    ch_bsj_bed_per_sample_tool_filtered = FILTER_BSJS( ch_bsj_bed_per_sample_tool, [] ).output
    ch_versions                         = ch_versions.mix(FILTER_BSJS.out.versions)

    //
    // MERGE BED FILES
    //

    BED_ADD_SAMPLE_TOOL( ch_bsj_bed_per_sample_tool_filtered, [] )
    ch_versions = ch_versions.mix(BED_ADD_SAMPLE_TOOL.out.versions)
    ch_bsj_bed_per_sample_tool_meta = BED_ADD_SAMPLE_TOOL.out.output

    COMBINE_TOOLS_PER_SAMPLE(
        ch_bsj_bed_per_sample_tool_meta
            .map{ meta, bed -> [ [id: meta.id], bed ] }
            .groupTuple(),
        params.max_shift,
        params.consider_strand,
        params.min_tools,
        1
    )
    ch_versions = ch_versions.mix(COMBINE_TOOLS_PER_SAMPLE.out.versions)
    ch_bsj_bed_per_sample = COMBINE_TOOLS_PER_SAMPLE.out.combined
        .filter{ meta, bed -> !bed.isEmpty() }

    COMBINE_SAMPLES(
        ch_bsj_bed_per_sample_tool_meta.map{ meta, bed -> [[id: "all"], bed] }.groupTuple(),
        params.max_shift,
        params.consider_strand,
        params.min_tools,
        params.min_samples
    )
    ch_versions = ch_versions.mix(COMBINE_SAMPLES.out.versions)
    ch_bsj_bed_combined = COMBINE_SAMPLES.out.combined
        .filter{ meta, bed -> !bed.isEmpty() }
        .collect()

    //
    // ANNOTATION
    //

    ADD_INTRONS(ch_gtf, [])
    ch_versions = ch_versions.mix(ADD_INTRONS.out.versions)

    EXTRACT_EXONS_INTRONS( ADD_INTRONS.out.gff, [] )
    ch_versions = ch_versions.mix(EXTRACT_EXONS_INTRONS.out.versions)

    ANNOTATE_COMBINED( ch_bsj_bed_combined, EXTRACT_EXONS_INTRONS.out.output, ch_annotation )
    ch_versions           = ch_versions.mix(ANNOTATE_COMBINED.out.versions)
    ch_bsj_bed12_combined = ANNOTATE_COMBINED.out.bed.collect()
    ch_bsj_gtf_combined   = ANNOTATE_COMBINED.out.gtf.collect()

    ANNOTATE_PER_SAMPLE( ch_bsj_bed_per_sample, EXTRACT_EXONS_INTRONS.out.output, ch_annotation )
    ch_versions             = ch_versions.mix(ANNOTATE_PER_SAMPLE.out.versions)
    ch_bsj_bed12_per_sample = ANNOTATE_PER_SAMPLE.out.bed
    ch_bsj_gtf_per_sample   = ANNOTATE_PER_SAMPLE.out.gtf

    ANNOTATE_PER_SAMPLE_TOOL( ch_bsj_bed_per_sample_tool, EXTRACT_EXONS_INTRONS.out.output, ch_annotation )
    ch_versions                  = ch_versions.mix(ANNOTATE_PER_SAMPLE_TOOL.out.versions)
    ch_bsj_bed12_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.bed
    ch_bsj_gtf_per_sample_tool   = ANNOTATE_PER_SAMPLE_TOOL.out.gtf

    //
    // FASTA WORKFLOW:
    //

    FASTA_COMBINED( ch_bsj_bed_combined, fasta )
    ch_versions = ch_versions.mix(FASTA_COMBINED.out.versions)
    ch_bsj_fasta_combined = FASTA_COMBINED.out.fasta

    FASTA_PER_SAMPLE( ch_bsj_bed_per_sample, fasta )
    ch_versions = ch_versions.mix(FASTA_PER_SAMPLE.out.versions)
    ch_bsj_fasta_per_sample = FASTA_PER_SAMPLE.out.fasta

    FASTA_PER_SAMPLE_TOOL( ch_bsj_bed_per_sample_tool, fasta )
    ch_versions = ch_versions.mix(FASTA_PER_SAMPLE_TOOL.out.versions)
    ch_bsj_fasta_per_sample_tool = FASTA_PER_SAMPLE_TOOL.out.fasta

    // STOP PIPELINE IF NO CIRCULAR RNAs WERE FOUND
    FAIL_ON_EMPTY(
        ch_bsj_bed_combined.ifEmpty([[id: "empty"], []]),
        // Make sure to wait for per-sample results
        Channel.empty()
            .mix(ch_bsj_bed12_combined)
            .mix(ch_bsj_bed12_per_sample)
            .mix(ch_bsj_bed12_per_sample_tool)
            .mix(ch_bsj_fasta_combined)
            .mix(ch_bsj_fasta_per_sample)
            .mix(ch_bsj_fasta_per_sample_tool)
            .map{ meta, f -> f }
            .collect()
    )

    emit:
    bed           = ch_bsj_bed_combined
    bed12         = ch_bsj_bed12_combined
    gtf           = ch_bsj_gtf_combined
    fasta         = ch_bsj_fasta_combined

    bed_per_sample_tool = ch_bsj_bed_per_sample_tool_meta

    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
