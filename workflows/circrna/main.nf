/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOWS:
include { paramsSummaryMap                 } from 'plugin/nf-validation'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { validateInputSamplesheet         } from '../../subworkflows/local/utils_nfcore_circrna_pipeline'

include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { PREPARE_GENOME                   } from '../../subworkflows/local/prepare_genome'
include { BSJ_DETECTION                    } from '../../subworkflows/local/bsj_detection'
include { ISOFORM_DETECTION                } from '../../subworkflows/local/isoform_detection'
include { ANNOTATION                       } from '../../subworkflows/local/annotation'
include { QUANTIFICATION                   } from '../../subworkflows/local/quantification'
include { MIRNA_PREDICTION                 } from '../../subworkflows/local/mirna_prediction'
include { STATISTICAL_TESTS                } from '../../subworkflows/local/statistical_tests'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES:
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
include { CAT_FASTQ                   } from '../../modules/nf-core/cat/fastq/main'

// SUBWORKFLOWS:
include { FASTQC_TRIMGALORE } from '../../subworkflows/nf-core/fastqc_trimgalore'
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
    ch_mature
    ch_annotation
    ch_versions
    ch_mirna

    main:

    ch_multiqc_files = Channel.empty()

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    ch_samplesheet
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            validateInputSamplesheet(it)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs ]
        }
        .set { ch_fastq }

    // MODULE:
    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (ch_fastq.multiple)
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // SUBORKFLOW:
    // Prepare index files &/or use iGenomes if chosen.
    PREPARE_GENOME (
        ch_fasta,
        ch_gtf
    )

    ch_gtf         = PREPARE_GENOME.out.gtf
    bowtie_index   = PREPARE_GENOME.out.bowtie
    bowtie2_index  = PREPARE_GENOME.out.bowtie2
    bwa_index      = PREPARE_GENOME.out.bwa
    chromosomes    = PREPARE_GENOME.out.chromosomes
    hisat2_index   = PREPARE_GENOME.out.hisat2
    star_index     = PREPARE_GENOME.out.star
    ch_versions    = ch_versions.mix(PREPARE_GENOME.out.versions)

    // MODULE: Run FastQC, trimgalore!
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files  = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    //
    // 2. BSJ Discovery
    //

    BSJ_DETECTION(
        FASTQC_TRIMGALORE.out.reads,
        ch_fasta,
        ch_gtf,
        ch_annotation,
        bowtie_index,
        bowtie2_index,
        bwa_index,
        chromosomes,
        hisat2_index,
        star_index,
        params.bsj_reads,
        params.exon_boundary
    )

    ch_multiqc_files  = ch_multiqc_files.mix(BSJ_DETECTION.out.multiqc_files)
    ch_versions = ch_versions.mix(BSJ_DETECTION.out.versions)

    //
    // 3. Isoform detection
    //

    ISOFORM_DETECTION(
        FASTQC_TRIMGALORE.out.reads,
        ch_cat_fastq,
        bwa_index,
        ch_fasta,
        ch_gtf,
        BSJ_DETECTION.out.bed12
    )

    //
    // 4. circRNA quantification
    //

    QUANTIFICATION(
        ch_gtf,
        ch_fasta,
        FASTQC_TRIMGALORE.out.reads,
        BSJ_DETECTION.out.bed12,
        BSJ_DETECTION.out.gtf,
        params.bootstrap_samples,
        ch_phenotype,
        PREPARE_GENOME.out.faidx
    )

    ch_versions = ch_versions.mix(QUANTIFICATION.out.versions)

    //
    // 5. miRNA prediction
    //

    if (params.mature) {
        MIRNA_PREDICTION(
            QUANTIFICATION.out.transcriptome,
            BSJ_DETECTION.out.bed12,
            ch_mature,
            ch_mirna,
            QUANTIFICATION.out.circular_tx_counts,
            QUANTIFICATION.out.rds
        )
        ch_versions = ch_versions.mix(MIRNA_PREDICTION.out.versions)
    }

    //
    // 6. Statistical tests
    //

    STATISTICAL_TESTS(
        QUANTIFICATION.out.se,
        QUANTIFICATION.out.gene_counts,
        QUANTIFICATION.out.circular_tx_counts,
        ch_phenotype
    )

    ch_versions = ch_versions.mix(STATISTICAL_TESTS.out.versions)


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    // MultiQC
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
