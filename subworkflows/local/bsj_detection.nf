// MODULES
include { GAWK as EXTRACT_COUNTS                             } from '../../modules/nf-core/gawk'
include { CSVTK_JOIN as COMBINE_COUNTS_PER_TOOL              } from '../../modules/nf-core/csvtk/join'
include { GAWK as FILTER_BSJS                                } from '../../modules/nf-core/gawk'
include { GAWK as BED_ADD_SAMPLE_TOOL                        } from '../../modules/nf-core/gawk'
include { COMBINEBEDS_READS                                  } from '../../modules/local/combinebeds/reads'
include { COMBINEBEDS_FILTER as COMBINE_TOOLS_PER_SAMPLE     } from '../../modules/local/combinebeds/filter'
include { COMBINEBEDS_SHIFTS as INVESTIGATE_SHIFTS           } from '../../modules/local/combinebeds/shifts'
include { COMBINEBEDS_FILTER as COMBINE_SAMPLES              } from '../../modules/local/combinebeds/filter'

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
    circexplorer2_index
    bsj_reads

    main:
    ch_versions                = Channel.empty()
    ch_bsj_bed_per_sample_tool = Channel.empty()
    ch_multiqc_files           = Channel.empty()
    fasta                      = ch_fasta.map{_meta, fasta -> fasta}
    gtf                        = ch_gtf.map{_meta, gtf -> gtf}

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
        CIRCEXPLORER2( fasta, STAR2PASS.out.junction, circexplorer2_index )
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
            STAR2PASS.out.junction, circexplorer2_index )
        ch_versions                = ch_versions.mix(MAPSPLICE.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(MAPSPLICE.out.bed)
    }

    ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool
        .filter{ _meta, bed -> !bed.isEmpty() }

    //
    // Analyze read-level agreement
    //

    tools_with_reads = ["find_circ", "segemehl", "dcc"]
    COMBINEBEDS_READS(
        ch_bsj_bed_per_sample_tool.filter{ _meta, _bed -> tools_with_reads.contains(_meta.tool) }
            .map{ meta, bed -> [[id: meta.id], meta.tool, bed] }
            .groupTuple()
    )
    ch_versions = ch_versions.mix(COMBINEBEDS_READS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COMBINEBEDS_READS.out.multiqc)

    //
    // QUANTIFY BSJs PER TOOL
    //

    EXTRACT_COUNTS( ch_bsj_bed_per_sample_tool, [], false )
    ch_versions = ch_versions.mix(EXTRACT_COUNTS.out.versions)

    COMBINE_COUNTS_PER_TOOL( EXTRACT_COUNTS.out.output
        .map{ meta, bed -> [[id: meta.tool], bed]}
        .groupTuple() )
    ch_versions = ch_versions.mix(COMBINE_COUNTS_PER_TOOL.out.versions)

    //
    // APPLY bsj_reads FILTER
    //

    ch_bsj_bed_per_sample_tool_filtered = FILTER_BSJS( ch_bsj_bed_per_sample_tool, [], false ).output
    ch_versions                         = ch_versions.mix(FILTER_BSJS.out.versions)

    //
    // MERGE BED FILES
    //

    BED_ADD_SAMPLE_TOOL( ch_bsj_bed_per_sample_tool_filtered, [], false )
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
        .filter{ _meta, bed -> !bed.isEmpty() }

    ch_all_samples = ch_bsj_bed_per_sample_tool_meta
        .map{ _meta, bed -> [[id: "all"], bed] }
        .groupTuple()

    COMBINE_SAMPLES(
        ch_all_samples,
        params.max_shift,
        params.consider_strand,
        params.min_tools,
        params.min_samples
    )
    ch_versions = ch_versions.mix(COMBINE_SAMPLES.out.versions)
    ch_bsj_bed_combined = COMBINE_SAMPLES.out.combined
        .filter{ _meta, bed -> !bed.isEmpty() }
        .collect()

    INVESTIGATE_SHIFTS(ch_all_samples)
    ch_versions = ch_versions.mix(INVESTIGATE_SHIFTS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(INVESTIGATE_SHIFTS.out.multiqc)

    //
    // ANNOTATION
    //

    ANNOTATE_COMBINED( ch_bsj_bed_combined, ch_annotation, fasta, circexplorer2_index )
    ch_versions           = ch_versions.mix(ANNOTATE_COMBINED.out.versions)
    ch_bsj_bed12_combined = ANNOTATE_COMBINED.out.bed12.collect()
    ch_bsj_gtf_combined   = ANNOTATE_COMBINED.out.gtf
    ch_bsj_fasta_combined = ANNOTATE_COMBINED.out.fasta

    ANNOTATE_PER_SAMPLE( ch_bsj_bed_per_sample, ch_annotation, fasta, circexplorer2_index )
    ch_versions             = ch_versions.mix(ANNOTATE_PER_SAMPLE.out.versions)
    ch_bsj_bed12_per_sample = ANNOTATE_PER_SAMPLE.out.bed12
    ch_bsj_fasta_per_sample = ANNOTATE_PER_SAMPLE.out.fasta

    ANNOTATE_PER_SAMPLE_TOOL( ch_bsj_bed_per_sample_tool, ch_annotation, fasta, circexplorer2_index )
    ch_versions                  = ch_versions.mix(ANNOTATE_PER_SAMPLE_TOOL.out.versions)
    ch_bsj_bed12_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.bed12
    ch_bsj_fasta_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.fasta

    // STOP PIPELINE IF NO CIRCULAR RNAs WERE FOUND
    Channel.empty()
        // First, make sure that all detection processes are finished
        // So that even if all hits are filtered out, we will still get the intermediate files
        .mix(ch_bsj_bed12_combined)
        .mix(ch_bsj_bed12_per_sample)
        .mix(ch_bsj_bed12_per_sample_tool)
        .mix(ch_bsj_fasta_combined)
        .mix(ch_bsj_fasta_per_sample)
        .mix(ch_bsj_fasta_per_sample_tool)
        // Then, check if any circular RNAs were found
        .map{ _meta, _f -> false }
        .mix( ch_bsj_bed_combined.map{ _meta, _f -> true } )
        // If no circular RNAs were found, stop the pipeline
        .filter{ it }
        .ifEmpty{
            error (
                "No circular RNAs were found by at least ${params.min_tools} tools and in at least ${params.min_samples} samples.\n" +
                "These thresholds can be adjusted using the parameters 'min_tools' and 'min_samples'.\n" +
                "Feel free to check the preliminary results in '${params.outdir}'\n" +
                (params.save_intermediates ? "" :
                "You can enable saving intermediate files by setting the parameter 'save_intermediates' to 'true'."))
        }

    emit:
    bed           = ch_bsj_bed_combined
    bed12         = ch_bsj_bed12_combined
    gtf           = ch_bsj_gtf_combined
    fasta         = ch_bsj_fasta_combined

    bed_per_sample_tool = ch_bsj_bed_per_sample_tool_meta

    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
