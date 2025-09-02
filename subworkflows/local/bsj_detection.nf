// MODULES
include { GAWK as EXTRACT_COUNTS                         } from '../../modules/nf-core/gawk'
include { CSVTK_JOIN as COMBINE_COUNTS_PER_TOOL          } from '../../modules/nf-core/csvtk/join'
include { GAWK as FILTER_BSJS                            } from '../../modules/nf-core/gawk'
include { GAWK as BED_ADD_SAMPLE_TOOL                    } from '../../modules/nf-core/gawk'
include { BEDTOOLS_INTERSECT as BLACKLIST                } from '../../modules/nf-core/bedtools/intersect'
include { COMBINEBEDS_READS                              } from '../../modules/local/combinebeds/reads'
include { COMBINEBEDS_FILTER as COMBINE_TOOLS_PER_SAMPLE } from '../../modules/local/combinebeds/filter'
include { COMBINEBEDS_SHIFTS as INVESTIGATE_SHIFTS       } from '../../modules/local/combinebeds/shifts'
include { COMBINEBEDS_FILTER as COMBINE_SAMPLES          } from '../../modules/local/combinebeds/filter'
include { FAIL_ON_EMPTY                                  } from '../../modules/local/fail_on_empty'

// SUBWORKFLOWS
include { SEGEMEHL                                       } from './detection_tools/segemehl'
include { STAR2PASS                                      } from './detection_tools/star2pass'
include { CIRCEXPLORER2                                  } from './detection_tools/circexplorer2'
include { CIRCRNA_FINDER                                 } from './detection_tools/circrna_finder'
include { FIND_CIRC                                      } from './detection_tools/find_circ'
include { CIRI                                           } from './detection_tools/ciri'
include { CIRCTOOLS                                      } from './detection_tools/circtools'
include { MAPSPLICE                                      } from './detection_tools/mapsplice'
include { PSIRC                                          } from './detection_tools/psirc'
include { ANNOTATION as ANNOTATE_COMBINED                } from './annotation'
include { ANNOTATION as ANNOTATE_PER_SAMPLE              } from './annotation'
include { ANNOTATION as ANNOTATE_PER_SAMPLE_TOOL         } from './annotation'

workflow BSJ_DETECTION {
    take:
    reads
    ch_fasta
    ch_gtf
    ch_blacklist
    ch_annotation
    bowtie_index
    bowtie2_index
    bwa_index
    chromosomes
    star_index
    circexplorer2_index
    psirc_index
    bsj_reads

    main:
    ch_versions = Channel.empty()
    ch_bsj_bed_per_sample_tool = Channel.empty()
    ch_multiqc_files = Channel.empty()
    fasta = ch_fasta.map { _meta, fasta -> fasta }
    gtf = ch_gtf.map { _meta, gtf -> gtf }

    def tools_selected = params.tools.split(',').collect { it.trim().toLowerCase() }
    def fli_tools_selected = params.fli_tools.split(',').collect { it.trim().toLowerCase() }

    // STAR 2-PASS-MODE
    star_ignore_sjdbgtf = true
    seq_center = params.seq_center ?: ''
    seq_platform = ''

    if ((tools_selected.intersect(['circexplorer2', 'circrna_finder', 'circtools', 'mapsplice']).size() > 0) || (fli_tools_selected.intersect(['circtools']).size() > 0)) {
        STAR2PASS(reads, star_index, ch_gtf, bsj_reads, star_ignore_sjdbgtf, seq_center, seq_platform)
        ch_versions = ch_versions.mix(STAR2PASS.out.versions)
        ch_star_bam = STAR2PASS.out.bam
    }
    else {
        ch_star_bam = Channel.empty()
    }

    //
    // DISCOVERY TOOLS:
    //

    if (tools_selected.size() == 0) {
        error('No tools selected for circRNA discovery.')
    }

    if (tools_selected.contains('segemehl')) {
        SEGEMEHL(reads, fasta, params.segemehl)
        ch_versions = ch_versions.mix(SEGEMEHL.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(SEGEMEHL.out.bed)
    }

    if (tools_selected.contains('circexplorer2')) {
        CIRCEXPLORER2(fasta, STAR2PASS.out.junction, circexplorer2_index)
        ch_versions = ch_versions.mix(CIRCEXPLORER2.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCEXPLORER2.out.bed)
    }

    if (tools_selected.contains('circrna_finder')) {
        CIRCRNA_FINDER(
            STAR2PASS.out.sam,
            STAR2PASS.out.junction,
            STAR2PASS.out.tab,
        )
        ch_versions = ch_versions.mix(CIRCRNA_FINDER.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCRNA_FINDER.out.bed)
    }

    if (tools_selected.contains('find_circ')) {
        FIND_CIRC(reads, bowtie2_index, ch_fasta)
        ch_versions = ch_versions.mix(FIND_CIRC.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(FIND_CIRC.out.bed)
    }

    if (tools_selected.contains('ciri')) {
        CIRI(reads, ch_fasta, ch_gtf, bwa_index)
        ch_versions = ch_versions.mix(CIRI.out.versions)
        ch_ciri_txt = CIRI.out.ciri_txt
        ch_ciri_sam = CIRI.out.ciri_sam
        ch_reads_fixed_length = CIRI.out.reads_fixed_length
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRI.out.bed)
    }
    else {
        ch_reads_fixed_length = Channel.empty()
        ch_ciri_txt = Channel.empty()
        ch_ciri_sam = Channel.empty()
    }

    if (tools_selected.contains('circtools')) {
        CIRCTOOLS(
            reads,
            ch_fasta,
            ch_gtf,
            star_index,
            STAR2PASS.out.junction,
            star_ignore_sjdbgtf,
            seq_platform,
            seq_center,
            bsj_reads,
        )
        ch_versions = ch_versions.mix(CIRCTOOLS.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(CIRCTOOLS.out.bed)
    }

    if (tools_selected.contains('mapsplice')) {
        MAPSPLICE(
            reads,
            gtf,
            fasta,
            bowtie_index,
            chromosomes,
            STAR2PASS.out.junction,
            circexplorer2_index,
        )
        ch_versions = ch_versions.mix(MAPSPLICE.out.versions)
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(MAPSPLICE.out.bed)
    }

    if (tools_selected.contains('psirc')) {
        PSIRC(reads, psirc_index)
        ch_versions = ch_versions.mix(PSIRC.out.versions)
        ch_psirc_bsj = PSIRC.out.output
        ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.mix(PSIRC.out.bed)
    }
    else {
        ch_psirc_bsj = Channel.empty()
    }

    ch_bsj_bed_per_sample_tool = ch_bsj_bed_per_sample_tool.filter { _meta, bed -> !bed.isEmpty() }

    if (params.blacklist) {
        BLACKLIST(ch_bsj_bed_per_sample_tool.combine(ch_blacklist), [[], []])
        ch_versions = ch_versions.mix(BLACKLIST.out.versions)
        ch_bsj_bed_per_sample_tool = BLACKLIST.out.intersect
    }

    //
    // Analyze read-level agreement
    //

    def tools_with_reads = ["find_circ", "segemehl", "circtools", "ciri"]
    def enabled_tools_with_reads = tools_selected.intersect(tools_with_reads)

    ch_bsj_bed_per_sample_tool_reads = ch_bsj_bed_per_sample_tool.filter { _meta, _bed -> tools_with_reads.contains(_meta.tool) }

    if (enabled_tools_with_reads.size() > 1) {
        COMBINEBEDS_READS(
            ch_bsj_bed_per_sample_tool_reads.map { meta, bed -> [[id: meta.id], meta.tool, bed] }.groupTuple()
        )
        ch_versions = ch_versions.mix(COMBINEBEDS_READS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(COMBINEBEDS_READS.out.multiqc)
        ch_bsj_reads = COMBINEBEDS_READS.out.combined
    }
    else {
        ch_bsj_reads = ch_bsj_bed_per_sample_tool_reads
    }

    //
    // QUANTIFY BSJs PER TOOL
    //

    EXTRACT_COUNTS(ch_bsj_bed_per_sample_tool, [], false)
    ch_versions = ch_versions.mix(EXTRACT_COUNTS.out.versions)

    COMBINE_COUNTS_PER_TOOL(
        EXTRACT_COUNTS.out.output.map { meta, bed -> [[id: meta.tool], bed] }.groupTuple()
    )
    ch_versions = ch_versions.mix(COMBINE_COUNTS_PER_TOOL.out.versions)

    //
    // APPLY bsj_reads FILTER
    //

    ch_bsj_bed_per_sample_tool_filtered = FILTER_BSJS(ch_bsj_bed_per_sample_tool, [], false).output
    ch_versions = ch_versions.mix(FILTER_BSJS.out.versions)

    //
    // MERGE BED FILES
    //

    BED_ADD_SAMPLE_TOOL(ch_bsj_bed_per_sample_tool_filtered, [], false)
    ch_versions = ch_versions.mix(BED_ADD_SAMPLE_TOOL.out.versions)
    ch_bsj_bed_per_sample_tool_meta = BED_ADD_SAMPLE_TOOL.out.output

    COMBINE_TOOLS_PER_SAMPLE(
        ch_bsj_bed_per_sample_tool_meta.map { meta, bed -> [[id: meta.id], bed] }.groupTuple(),
        params.max_shift,
        params.consider_strand,
        params.min_tools,
        1,
    )
    ch_versions = ch_versions.mix(COMBINE_TOOLS_PER_SAMPLE.out.versions)
    ch_bsj_bed_per_sample = COMBINE_TOOLS_PER_SAMPLE.out.combined.filter { _meta, bed -> !bed.isEmpty() }

    ch_all_samples = ch_bsj_bed_per_sample_tool_meta
        .map { _meta, bed -> [[id: "all"], bed] }
        .groupTuple()

    COMBINE_SAMPLES(
        ch_all_samples,
        params.max_shift,
        params.consider_strand,
        params.min_tools,
        params.min_samples,
    )
    ch_versions = ch_versions.mix(COMBINE_SAMPLES.out.versions)
    ch_bsj_bed_combined = COMBINE_SAMPLES.out.combined
        .filter { _meta, bed -> !bed.isEmpty() }
        .collect()

    INVESTIGATE_SHIFTS(ch_all_samples)
    ch_versions = ch_versions.mix(INVESTIGATE_SHIFTS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(INVESTIGATE_SHIFTS.out.multiqc)

    //
    // ANNOTATION
    //

    ANNOTATE_COMBINED(ch_bsj_bed_combined, ch_annotation, fasta, circexplorer2_index)
    ch_versions = ch_versions.mix(ANNOTATE_COMBINED.out.versions)
    ch_bsj_bed12_combined = ANNOTATE_COMBINED.out.bed12.collect()
    ch_bsj_gtf_combined = ANNOTATE_COMBINED.out.gtf
    ch_bsj_fasta_combined = ANNOTATE_COMBINED.out.fasta

    ANNOTATE_PER_SAMPLE(ch_bsj_bed_per_sample, ch_annotation, fasta, circexplorer2_index)
    ch_versions = ch_versions.mix(ANNOTATE_PER_SAMPLE.out.versions)
    ch_bsj_bed12_per_sample = ANNOTATE_PER_SAMPLE.out.bed12
    ch_bsj_fasta_per_sample = ANNOTATE_PER_SAMPLE.out.fasta

    ANNOTATE_PER_SAMPLE_TOOL(ch_bsj_bed_per_sample_tool, ch_annotation, fasta, circexplorer2_index)
    ch_versions = ch_versions.mix(ANNOTATE_PER_SAMPLE_TOOL.out.versions)
    ch_bsj_bed12_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.bed12
    ch_bsj_fasta_per_sample_tool = ANNOTATE_PER_SAMPLE_TOOL.out.fasta

    // STOP PIPELINE IF NO CIRCULAR RNAs WERE FOUND
    FAIL_ON_EMPTY(
        ch_bsj_bed_combined.ifEmpty([[id: "empty"], []]),
        Channel.empty().mix(ch_bsj_bed12_combined).mix(ch_bsj_bed12_per_sample).mix(ch_bsj_bed12_per_sample_tool).mix(ch_bsj_fasta_combined).mix(ch_bsj_fasta_per_sample).mix(ch_bsj_fasta_per_sample_tool).map { _meta, f -> f }.collect(),
    )

    emit:
    bed                 = ch_bsj_bed_combined
    bed12               = ch_bsj_bed12_combined
    gtf                 = ch_bsj_gtf_combined
    fasta               = ch_bsj_fasta_combined
    bed_reads           = ch_bsj_reads
    bed_per_sample      = ch_bsj_bed_per_sample
    bed_per_sample_tool = ch_bsj_bed_per_sample_tool_meta
    ciri_txt            = ch_ciri_txt
    ciri_sam            = ch_ciri_sam
    reads_fixed_length  = ch_reads_fixed_length
    psirc_bsj           = ch_psirc_bsj
    star_bam            = ch_star_bam
    star_junction       = STAR2PASS.out.junction
    multiqc_files       = ch_multiqc_files
    versions            = ch_versions
}
