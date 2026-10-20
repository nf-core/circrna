// MODULES
include { GAWK as EXTRACT_COUNTS                         } from '../../../modules/nf-core/gawk'
include { GAWK as FILTER_BSJS                            } from '../../../modules/nf-core/gawk'
include { GAWK as BED_ADD_SAMPLE_TOOL                    } from '../../../modules/nf-core/gawk'
include { BEDTOOLS_INTERSECT as BLACKLIST                } from '../../../modules/nf-core/bedtools/intersect'
include { JOIN_SAMPLES as COMBINE_COUNTS_PER_TOOL        } from '../../../modules/local/matrix/join_samples'
include { COMBINEBEDS_READS                              } from '../../../modules/local/combinebeds/reads'
include { COMBINEBEDS_FILTER as COMBINE_TOOLS_PER_SAMPLE } from '../../../modules/local/combinebeds/filter'
include { COMBINEBEDS_SHIFTS as INVESTIGATE_SHIFTS       } from '../../../modules/local/combinebeds/shifts'
include { COMBINEBEDS_FILTER as COMBINE_SAMPLES          } from '../../../modules/local/combinebeds/filter'
include { FAIL_ON_EMPTY                                  } from '../../../modules/local/fail_on_empty'

// SUBWORKFLOWS
include { SEGEMEHL                                       } from '../detection_tools/segemehl'
include { STAR2PASS                                      } from '../detection_tools/star2pass'
include { CIRCEXPLORER2                                  } from '../detection_tools/circexplorer2'
include { CIRCRNA_FINDER                                 } from '../detection_tools/circrna_finder'
include { FIND_CIRC                                      } from '../detection_tools/find_circ'
include { CIRI2                                          } from '../detection_tools/ciri2'
include { CIRI3                                          } from '../detection_tools/ciri3'
include { CIRCTOOLS                                      } from '../detection_tools/circtools'
include { MAPSPLICE                                      } from '../detection_tools/mapsplice'
include { PSIRC                                          } from '../detection_tools/psirc'
include { CIRCMINER                                      } from '../detection_tools/circminer'
include { ANNOTATION as ANNOTATE_COMBINED                } from '../annotation'
include { ANNOTATION as ANNOTATE_PER_SAMPLE              } from '../annotation'
include { ANNOTATION as ANNOTATE_PER_SAMPLE_TOOL         } from '../annotation'
include { starBasedTools                                 } from '../utils_nfcore_circrna_pipeline'
include { readCarryingTools                              } from '../utils_nfcore_circrna_pipeline'

workflow BSJ_DETECTION {
    take:
    reads
    ch_fasta
    ch_gtf
    blacklist
    ch_annotation
    bowtie_index
    bowtie2_index
    ch_bwa_sam
    chromosomes
    star_index
    circexplorer2_index
    psirc_index
    bsj_reads
    tools
    fli_tools
    segemehl_index
    max_shift
    consider_strand
    min_tools
    min_samples
    exons_only
    ch_faidx
    circminer_index
    circminer_gtf
    ch_reads_fixed_length

    main:
    ch_bsj_bed_per_sample_tool = channel.empty()
    ch_multiqc_files = channel.empty()
    fasta = ch_fasta.map { _meta, fasta -> fasta }
    gtf = ch_gtf.map { _meta, gtf -> gtf }

    // STAR 2-PASS-MODE
    star_ignore_sjdbgtf = true

    if (starBasedTools(tools).size() > 0 || fli_tools.contains('circtools')) {
        STAR2PASS(reads, star_index, ch_gtf, bsj_reads, star_ignore_sjdbgtf)
        ch_star_bam = STAR2PASS.out.bam
        ch_star_junction = STAR2PASS.out.junction
    }
    else {
        ch_star_bam = channel.empty()
        ch_star_junction = channel.empty()
    }

    //
    // DISCOVERY TOOLS:
    //

    if (tools.size() == 0) {
        error('No tools selected for circRNA discovery.')
    }

    if (tools.contains('segemehl')) {
        ch_segemehl_index = segemehl_index
            ? channel.value([[id: 'segemehl_reference'], file(segemehl_index)])
            : null
        SEGEMEHL(reads, ch_fasta, ch_segemehl_index)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(SEGEMEHL.out.bed)
    }

    if (tools.contains('circexplorer2')) {
        CIRCEXPLORER2(fasta, STAR2PASS.out.junction, circexplorer2_index)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCEXPLORER2.out.bed)
    }

    if (tools.contains('circrna_finder')) {
        CIRCRNA_FINDER(
            STAR2PASS.out.sam,
            STAR2PASS.out.junction,
            STAR2PASS.out.tab,
        )
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCRNA_FINDER.out.bed)
    }

    if (tools.contains('find_circ')) {
        FIND_CIRC(reads, bowtie2_index, ch_fasta)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(FIND_CIRC.out.bed)
    }

    if (tools.contains('ciri2')) {
        CIRI2(ch_bwa_sam, ch_fasta, ch_gtf)
        ch_ciri2_txt = CIRI2.out.ciri2_txt
        ch_ciri2_sam = CIRI2.out.ciri2_sam
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRI2.out.bed)
    }
    else {
        ch_ciri2_txt = channel.empty()
        ch_ciri2_sam = channel.empty()
    }

    if (tools.contains('ciri3')) {
        CIRI3(ch_bwa_sam, ch_fasta, ch_gtf)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRI3.out.bed)
    }

    if (tools.contains('circtools')) {
        CIRCTOOLS(
            reads,
            ch_fasta,
            ch_gtf,
            star_index,
            STAR2PASS.out.junction,
            star_ignore_sjdbgtf,
            bsj_reads,
        )
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCTOOLS.out.bed)
    }

    if (tools.contains('mapsplice')) {
        MAPSPLICE(
            reads,
            gtf,
            fasta,
            bowtie_index,
            chromosomes,
            STAR2PASS.out.junction,
            circexplorer2_index,
        )
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(MAPSPLICE.out.bed)
    }

    if (tools.contains('psirc')) {
        PSIRC(reads, psirc_index, ch_fasta, ch_faidx)
        ch_psirc_bsj = PSIRC.out.output
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(PSIRC.out.bed)
    }
    else {
        ch_psirc_bsj = channel.empty()
    }

    if (tools.contains('circminer')) {
        CIRCMINER(reads, circminer_index, circminer_gtf, ch_fasta, ch_faidx)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCMINER.out.bed)
    }

    ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.filter { _meta, bed -> !bed.isEmpty() }

    if (blacklist) {
        BLACKLIST(ch_bsj_bed_per_sample_tool.combine(channel.value(blacklist)), [[], []], [])
        ch_bsj_bed_per_sample_tool = BLACKLIST.out.intersect
    }

    //
    // Analyze read-level agreement
    //

    def enabled_tools_with_reads = readCarryingTools(tools)

    ch_bsj_bed_per_sample_tool_reads = ch_bsj_bed_per_sample_tool.filter { _meta, _bed -> readCarryingTools([_meta.tool]).size() > 0 }

    if (enabled_tools_with_reads.size() > 1) {
        COMBINEBEDS_READS(
            ch_bsj_bed_per_sample_tool_reads.map { meta, bed -> [[id: meta.id], meta.tool, bed] }.groupTuple(),
            consider_strand,
            2,
        )
        ch_multiqc_files = ch_multiqc_files.mix(COMBINEBEDS_READS.out.multiqc)
        ch_bsj_reads = COMBINEBEDS_READS.out.combined
    }
    else {
        ch_bsj_reads = ch_bsj_bed_per_sample_tool_reads.map { meta, bed -> [[id: meta.id], bed] }
    }

    //
    // QUANTIFY BSJs PER TOOL
    //

    EXTRACT_COUNTS(ch_bsj_bed_per_sample_tool.map { meta, bed -> [meta, bed, "counts.tsv"] }, [], false)

    COMBINE_COUNTS_PER_TOOL(
        EXTRACT_COUNTS.out.output.map { meta, counts -> [[id: meta.tool], meta.id, counts] }.groupTuple(),
        "id",
    )

    //
    // APPLY bsj_reads FILTER
    //

    ch_bsj_bed_per_sample_tool_filtered = FILTER_BSJS(ch_bsj_bed_per_sample_tool.map { meta, bed -> [meta, bed, "${meta.tool}.filtered.bed"] }, [], false).output

    //
    // MERGE BED FILES
    //

    BED_ADD_SAMPLE_TOOL(ch_bsj_bed_per_sample_tool_filtered.map { meta, bed -> [meta, bed, "meta.bed"] }, [], false)
    ch_bsj_bed_per_sample_tool_meta = BED_ADD_SAMPLE_TOOL.out.output

    COMBINE_TOOLS_PER_SAMPLE(
        ch_bsj_bed_per_sample_tool_meta.map { meta, bed -> [[id: meta.id], bed] }.groupTuple(),
        max_shift,
        consider_strand,
        min_tools,
        1,
        "combined.bed",
    )
    ch_bsj_bed_per_sample = COMBINE_TOOLS_PER_SAMPLE.out.combined.filter { _meta, bed -> !bed.isEmpty() }

    ch_all_samples = ch_bsj_bed_per_sample_tool_meta
        .map { _meta, bed -> [[id: "all"], bed] }
        .groupTuple()

    COMBINE_SAMPLES(
        ch_all_samples,
        max_shift,
        consider_strand,
        min_tools,
        min_samples,
        "combined.bed",
    )
    ch_bsj_bed_combined = COMBINE_SAMPLES.out.combined
        .filter { _meta, bed -> !bed.isEmpty() }
        .collect()

    INVESTIGATE_SHIFTS(ch_all_samples)
    ch_multiqc_files = ch_multiqc_files.mix(INVESTIGATE_SHIFTS.out.multiqc.map { _meta, multiqc -> multiqc })

    //
    // ANNOTATION
    //

    ANNOTATE_COMBINED(ch_bsj_bed_combined, ch_annotation, fasta, circexplorer2_index, exons_only)
    ch_bsj_bed12_combined = ANNOTATE_COMBINED.out.bed12.collect()
    ch_bsj_gtf_combined = ANNOTATE_COMBINED.out.gtf
    ch_bsj_fasta_combined = ANNOTATE_COMBINED.out.fasta

    ANNOTATE_PER_SAMPLE(ch_bsj_bed_per_sample, ch_annotation, fasta, circexplorer2_index, exons_only)
    ch_bsj_bed12_per_sample = ANNOTATE_PER_SAMPLE.out.bed12
    ch_bsj_fasta_per_sample = ANNOTATE_PER_SAMPLE.out.fasta

    ANNOTATE_PER_SAMPLE_TOOL(ch_bsj_bed_per_sample_tool, ch_annotation, fasta, circexplorer2_index, exons_only)
    ch_bsj_bed12_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.bed12
    ch_bsj_fasta_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.fasta

    // STOP PIPELINE IF NO CIRCULAR RNAs WERE FOUND
    def fail_message = "No circular RNAs were found by at least ${min_tools} tools and in at least ${min_samples} samples.\n" + "Feel free to check the preliminary results in the output directory.\n" + "If they were not saved, you can enable saving intermediate files by setting the parameter 'save_intermediates' to 'true'."

    FAIL_ON_EMPTY(
        ch_bsj_bed_combined.ifEmpty([[id: "empty"], []]),
        channel.empty().mix(ch_bsj_bed12_combined).mix(ch_bsj_bed12_per_sample).mix(ch_bsj_bed12_per_sample_tool).mix(ch_bsj_fasta_combined).mix(ch_bsj_fasta_per_sample).mix(ch_bsj_fasta_per_sample_tool).map { _meta, f -> f }.collect().ifEmpty([]),
        fail_message,
    )

    emit:
    bed                 = ch_bsj_bed_combined
    bed12               = ch_bsj_bed12_combined
    gtf                 = ch_bsj_gtf_combined
    fasta               = ch_bsj_fasta_combined
    bed_reads           = ch_bsj_reads
    bed_per_sample      = ch_bsj_bed_per_sample
    bed12_per_sample    = ch_bsj_bed12_per_sample
    bed_per_sample_tool = ch_bsj_bed_per_sample_tool_meta
    ciri2_txt           = ch_ciri2_txt
    ciri2_sam           = ch_ciri2_sam
    reads_fixed_length  = ch_reads_fixed_length
    psirc_bsj           = ch_psirc_bsj
    star_bam            = ch_star_bam
    star_junction       = ch_star_junction
    multiqc_files       = ch_multiqc_files.flatten()
}
