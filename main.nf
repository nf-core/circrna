#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/circrna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/circrna
    Website: https://nf-co.re/circrna
    Slack  : https://nfcore.slack.com/channels/circrna
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CIRCRNA                 } from './workflows/circrna.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_circrna_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_circrna_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_circrna_pipeline'
include { parseToolList           } from './subworkflows/local/utils_nfcore_circrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta     = getGenomeAttribute('fasta')
params.gtf       = getGenomeAttribute('gtf')
params.bwa       = getGenomeAttribute('bwa')
params.star      = getGenomeAttribute('star')
params.bowtie    = getGenomeAttribute('bowtie')
params.bowtie2   = getGenomeAttribute('bowtie2')
params.mature    = getGenomeAttribute('mature')
params.blacklist = getGenomeAttribute('blacklist')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CIRCRNA {
    take:
    ch_samplesheet

    main:


    //
    // WORKFLOW: Run nf-core/circrna workflow
    //
    ch_fasta = channel.value([[id: "fasta"], file(params.fasta, checkIfExists: true)])
    ch_gtf = channel.value([[id: "gtf"], file(params.gtf, checkIfExists: true)])
    blacklist = params.blacklist ? file(params.blacklist, checkIfExists: true) : null
    mature = params.mature ? [[id: "mature"], file(params.mature, checkIfExists: true)] : null
    mirna_expression = params.mature && params.mirna_expression ? [[id: "mirna"], file(params.mirna_expression, checkIfExists: true)] : null
    ch_phenotype = params.phenotype ? channel.value([[id: "phenotype"], file(params.phenotype, checkIfExists: true)]) : channel.empty()
    ch_annotation = params.annotation && params.annotation instanceof String
        ? channel.fromSamplesheet("annotation", parameters_schema: "${projectDir}/nextflow_schema.json")
        : channel.empty()

    CIRCRNA(
        ch_samplesheet,
        ch_phenotype,
        ch_fasta,
        ch_gtf,
        blacklist,
        mature,
        ch_annotation,
        mirna_expression,
        parseToolList(params.tools),
        parseToolList(params.fli_tools),
        parseToolList(params.quantification_tools),
        parseToolList(params.mirna_tools),
        params.bowtie ? file(params.bowtie, checkIfExists: true) : null,
        params.bowtie2 ? file(params.bowtie2, checkIfExists: true) : null,
        params.bwa ? file(params.bwa, checkIfExists: true) : null,
        params.hisat2 ? file(params.hisat2, checkIfExists: true) : null,
        params.star ? file(params.star, checkIfExists: true) : null,
        params.segemehl ? file(params.segemehl, checkIfExists: true) : null,
        params.circminer ? file(params.circminer, checkIfExists: true) : null,
        params.longread,
        params.skip_fastqc,
        params.skip_trimming,
        params.bsj_reads,
        params.max_shift,
        params.consider_strand,
        params.min_tools,
        params.min_samples,
        params.exons_only,
        params.bootstrap_samples,
        params.mirna_min_tools,
        params.mirna_min_reads,
        params.mirna_min_sample_percentage,
        params.multiqc_config ? file(params.multiqc_config) : null,
        params.multiqc_logo ? file(params.multiqc_logo) : null,
        params.outdir,
    )

    emit:
    multiqc_report = CIRCRNA.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CIRCRNA(
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_CIRCRNA.out.multiqc_report,
    )
}
