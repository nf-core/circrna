/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome'
include { BSJ_DETECTION                    } from '../subworkflows/local/bsj_detection'
include { FLI_DETECTION                    } from '../subworkflows/local/fli_detection'
include { COMBINE_TRANSCRIPTOMES           } from '../subworkflows/local/combine_transcriptomes'
include { QUANTIFICATION                   } from '../subworkflows/local/quantification'
include { MIRNA_PREDICTION                 } from '../subworkflows/local/mirna_prediction'
include { STATISTICAL_TESTS                } from '../subworkflows/local/statistical_tests'
include { LONGREAD                         } from '../subworkflows/local/longread'
include { CIRI_READS                       } from '../subworkflows/local/ciri_reads'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES:
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { BWA_MEM as CIRI_BWA_MEM     } from '../modules/nf-core/bwa/mem/main'

// SUBWORKFLOWS:
include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CIRCRNA {
    take:
    ch_samplesheet
    ch_phenotype
    ch_fasta
    ch_gtf
    blacklist                   // path or null
    mature                      // [meta, path] or null
    ch_annotation
    ch_versions
    mirna_expression            // [meta, path] or null

    // tool selection
    tools                       // List<String>
    fli_tools                   // List<String>
    quantification_tools        // List<String>
    mirna_tools                 // List<String>

    // prebuilt indices (path or null)
    bowtie
    bowtie2
    bwa
    hisat2
    star
    segemehl
    circminer

    // settings
    longread
    skip_fastqc
    skip_trimming
    bsj_reads
    max_shift
    consider_strand
    min_tools
    min_samples
    exons_only
    bootstrap_samples
    mirna_min_tools
    mirna_min_reads
    mirna_min_sample_percentage

    // reporting
    multiqc_config              // path or null
    multiqc_logo                // path or null
    outdir

    main:

    ch_multiqc_files = channel.empty()

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    ch_samplesheet
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    // MODULE:
    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (ch_fastq.multiple)
        .reads
        .mix(
            ch_fastq.single
        )
        .set { ch_cat_fastq }

    // SUBORKFLOW:
    // Prepare index files &/or use iGenomes if chosen.
    PREPARE_GENOME (
        ch_fasta,
        ch_gtf,
        tools,
        quantification_tools,
        bowtie,
        bowtie2,
        bwa,
        hisat2,
        star,
        fli_tools,
        circminer
    )

    ch_gtf              = PREPARE_GENOME.out.gtf
    bowtie_index        = PREPARE_GENOME.out.bowtie
    bowtie2_index       = PREPARE_GENOME.out.bowtie2
    bwa_index           = PREPARE_GENOME.out.bwa
    chromosomes         = PREPARE_GENOME.out.chromosomes
    circexplorer2_index = PREPARE_GENOME.out.circexplorer2
    star_index          = PREPARE_GENOME.out.star
    psirc_index         = PREPARE_GENOME.out.psirc
    circminer_index     = PREPARE_GENOME.out.circminer
    ch_versions         = ch_versions.mix(PREPARE_GENOME.out.versions)



    if (longread) {
        LONGREAD(
            ch_cat_fastq,
            "hg38"
        )
        ch_versions = ch_versions.mix(LONGREAD.out.versions)
    } else {
        // MODULE: Run FastQC, trimgalore!
        FASTQC_TRIMGALORE (
            ch_cat_fastq,
            skip_fastqc,
            skip_trimming
        )
        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
        ch_multiqc_files  = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{ _meta, zip -> zip }.ifEmpty([]))
        ch_multiqc_files  = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_log.collect{ _meta, log -> log }.ifEmpty([]))

        //
        // 2. BSJ Discovery
        //

        if (tools.intersect(['ciri2', 'ciri3']).size() > 0 || quantification_tools.contains('ciri3')) {
            if (fli_tools.contains('cirifull')) {
                CIRI_READS(FASTQC_TRIMGALORE.out.reads)
                ch_bwa_reads = CIRI_READS.out.reads
                ch_reads_fixed_length = CIRI_READS.out.reads
                ch_versions = ch_versions.mix(CIRI_READS.out.versions)
            }
            else {
                ch_bwa_reads = FASTQC_TRIMGALORE.out.reads
                ch_reads_fixed_length = channel.empty()
            }

            CIRI_BWA_MEM(ch_bwa_reads, bwa_index.first(), ch_fasta, true)
            ch_bwa_sam = CIRI_BWA_MEM.out.sam
        }
        else {
            ch_bwa_sam = channel.empty()
            ch_reads_fixed_length = channel.empty()
        }

        BSJ_DETECTION(
            FASTQC_TRIMGALORE.out.reads,
            ch_fasta,
            ch_gtf,
            blacklist,
            ch_annotation,
            bowtie_index,
            bowtie2_index,
            ch_bwa_sam,
            chromosomes,
            star_index,
            circexplorer2_index,
            psirc_index,
            bsj_reads,
            tools,
            fli_tools,
            segemehl,
            max_shift,
            consider_strand,
            min_tools,
            min_samples,
            exons_only,
            PREPARE_GENOME.out.faidx,
            circminer_index,
            PREPARE_GENOME.out.circminer_gtf,
            ch_reads_fixed_length
        )

        ch_multiqc_files  = ch_multiqc_files.mix(BSJ_DETECTION.out.multiqc_files)
        ch_versions = ch_versions.mix(BSJ_DETECTION.out.versions)

        //
        // 3. FLI Detection
        //

        FLI_DETECTION(
            FASTQC_TRIMGALORE.out.reads,
            BSJ_DETECTION.out.reads_fixed_length,
            ch_fasta,
            ch_gtf,
            bwa_index,
            BSJ_DETECTION.out.ciri2_txt,
            BSJ_DETECTION.out.ciri2_sam,
            BSJ_DETECTION.out.bed12,
            BSJ_DETECTION.out.bed12_per_sample,
            BSJ_DETECTION.out.bed_reads,
            psirc_index,
            BSJ_DETECTION.out.psirc_bsj,
            BSJ_DETECTION.out.star_bam,
            BSJ_DETECTION.out.star_junction,
            BSJ_DETECTION.out.bed_per_sample,
            fli_tools,
            max_shift
        )
        ch_versions = ch_versions.mix(FLI_DETECTION.out.versions)

        COMBINE_TRANSCRIPTOMES(
            ch_fasta,
            ch_gtf,
            BSJ_DETECTION.out.gtf
        )

        ch_versions = ch_versions.mix(COMBINE_TRANSCRIPTOMES.out.versions)

        //
        // 4. circRNA quantification
        //

        QUANTIFICATION(
            FASTQC_TRIMGALORE.out.reads,
            ch_gtf,
            ch_fasta,
            COMBINE_TRANSCRIPTOMES.out.fasta,
            COMBINE_TRANSCRIPTOMES.out.gtf,
            BSJ_DETECTION.out.bed12,
            BSJ_DETECTION.out.gtf,
            BSJ_DETECTION.out.bed_per_sample_tool,
            bootstrap_samples,
            ch_phenotype,
            PREPARE_GENOME.out.faidx,
            PREPARE_GENOME.out.bwa,
            PREPARE_GENOME.out.hisat2,
            quantification_tools,
            max_shift,
            consider_strand,
            ch_bwa_sam
        )

        ch_versions = ch_versions.mix(QUANTIFICATION.out.versions)

        //
        // 5. miRNA prediction
        //

        if (mature) {
            MIRNA_PREDICTION(
                COMBINE_TRANSCRIPTOMES.out.fasta,
                BSJ_DETECTION.out.bed12,
                channel.value(mature),
                mirna_expression,
                QUANTIFICATION.out.circ,
                QUANTIFICATION.out.rds,
                mirna_tools,
                mirna_min_tools,
                mirna_min_reads,
                mirna_min_sample_percentage
            )
            ch_versions = ch_versions.mix(MIRNA_PREDICTION.out.versions)
        }

        //
        // 6. Statistical tests
        //

        STATISTICAL_TESTS(
            QUANTIFICATION.out.gene,
            QUANTIFICATION.out.circ,
            QUANTIFICATION.out.ciriquant,
            QUANTIFICATION.out.stringtie,
            ch_phenotype
        )
        ch_versions = ch_versions.mix(STATISTICAL_TESTS.out.versions)

    }



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'nf_core_circrna_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    // MultiQC
    ch_multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = multiqc_config ? channel.fromPath(multiqc_config) : channel.empty()
    ch_multiqc_logo          = multiqc_logo   ? channel.fromPath(multiqc_logo)   : channel.empty()
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    ch_multiqc_config_files = ch_multiqc_config
        .mix(ch_multiqc_custom_config)
        .collect()
        .map { configs -> [configs] }   // double-wrap so combine appends [configs] as one element

    ch_multiqc_logo_files = ch_multiqc_logo
        .collect()
        .map { logos -> [logos] }       // double-wrap so combine appends [logos] as one element

    MULTIQC (
        ch_multiqc_files
            .collect()
            .map { files -> [[id: "multiqc"], files] }
            .combine(ch_multiqc_config_files)
            .combine(ch_multiqc_logo_files)
            .map { meta, files, configs, logos ->
                [meta, files, configs, logos, [], []]
            }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
