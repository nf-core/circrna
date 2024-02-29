/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_phenotype   = params.phenotype && params.module.contains('differential_expression') ? file(params.phenotype, checkIfExists:true) : Channel.empty()
ch_fasta       = params.fasta ? file(params.fasta) : 'null'
ch_gtf         = params.gtf ? file(params.gtf) : 'null'
ch_mature      = params.mature && params.module.contains('mirna_prediction') ? file(params.mature) : Channel.empty()

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
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { PREPARE_GENOME                   } from '../../subworkflows/local/prepare_genome'
include { CIRCRNA_DISCOVERY                } from '../../subworkflows/local/circrna_discovery'
include { MIRNA_PREDICTION                 } from '../../subworkflows/local/mirna_prediction'
include { DIFFERENTIAL_EXPRESSION          } from '../../subworkflows/local/differential_expression'

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
    ch_versions

    main:

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    ch_samplesheet
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

    // MODULE:
    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // SUBORKFLOW:
    // Prepare index files &/or use iGenomes if chosen.
    PREPARE_GENOME (
        ch_fasta,
        ch_gtf
    )

    // Stage the indices via newly created indices, iGenomes or empty list if tool not selected.
    bowtie_index   = params.fasta ? params.bowtie ? Channel.fromPath(params.bowtie) : PREPARE_GENOME.out.bowtie : []
    bowtie2_index  = params.fasta ? params.bowtie2 ? Channel.fromPath(params.bowtie2).map{ it -> [[id:'bowtie2'], it] } : PREPARE_GENOME.out.bowtie2 : []
    bwa_index      = params.fasta ? params.bwa ? Channel.fromPath(params.bwa).map{ it -> [[id:'bwa'], it] } : PREPARE_GENOME.out.bwa : []
    chromosomes    = params.fasta && ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) ? PREPARE_GENOME.out.chromosomes : []
    hisat2_index   = params.fasta ? params.hisat2 && ( params.tool.contains('ciriquant') || params.module.contains('differential_expression') ) ? Channel.fromPath(params.hisat2).map{ [[id:"hisat2"], it]} : PREPARE_GENOME.out.hisat2 : []
    star_index     = params.fasta ? params.star ? Channel.fromPath(params.star).map{[[id:'star'], it]}: PREPARE_GENOME.out.star : []
    segemehl_index = params.fasta ? params.segemehl ? Channel.fromPath(params.segemehl) : PREPARE_GENOME.out.segemehl : []
    ch_versions    = ch_versions.mix(PREPARE_GENOME.out.versions)

    // MODULE: Run FastQC, trimgalore!
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    //
    // 2. circRNA Discovery
    //

    CIRCRNA_DISCOVERY(
        FASTQC_TRIMGALORE.out.reads,
        ch_fasta,
        ch_gtf,
        bowtie_index,
        bowtie2_index,
        bwa_index,
        chromosomes,
        hisat2_index,
        segemehl_index,
        star_index,
        params.bsj_reads,
        params.tool_filter,
        params.duplicates_fun,
        params.exon_boundary
    )

    ch_versions = ch_versions.mix(CIRCRNA_DISCOVERY.out.versions)

    //
    // 3. miRNA prediction
    //

    MIRNA_PREDICTION(
        CIRCRNA_DISCOVERY.out.fasta,
        CIRCRNA_DISCOVERY.out.circrna_bed12,
        ch_mature
    )

    ch_versions = ch_versions.mix(MIRNA_PREDICTION.out.versions)

    //
    // 4. Differential expression tests
    //

    ch_ensembl_database_map = params.module.contains('differential_expression') ? Channel.fromPath("${projectDir}/bin/ensembl_database_map.txt") : Channel.empty()

    DIFFERENTIAL_EXPRESSION(
        FASTQC_TRIMGALORE.out.reads,
        ch_gtf,
        ch_fasta,
        hisat2_index,
        PREPARE_GENOME.out.splice_sites,
        ch_phenotype,
        CIRCRNA_DISCOVERY.out.dea_matrix,
        CIRCRNA_DISCOVERY.out.clr_matrix,
        params.species,
        ch_ensembl_database_map,
        params.exon_boundary
    )

    ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)
    ch_reports  = ch_reports.mix(DIFFERENTIAL_EXPRESSION.out.reports)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    // MODULE: MultiQC
    workflow_summary    = WorkflowCircrna.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCircrna.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))

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
