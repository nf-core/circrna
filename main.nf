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

include { CIRCRNA                 } from './workflows/circrna'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_circrna_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_circrna_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_circrna_pipeline'

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

    ch_versions = channel.empty()

    //
    // WORKFLOW: Run nf-core/circrna workflow
    //
    ch_fasta = channel.value([[id: "fasta"], file(params.fasta, checkIfExists: true)])
    ch_gtf = channel.value([[id: "gtf"], file(params.gtf, checkIfExists: true)])
    ch_blacklist = params.blacklist ? channel.value(file(params.blacklist, checkIfExists: true)) : channel.empty()
    ch_mature = params.mature ? channel.value([[id: "mature"], file(params.mature, checkIfExists: true)]) : channel.empty()
    ch_phenotype = params.phenotype ? channel.value([[id: "phenotype"], file(params.phenotype, checkIfExists: true)]) : channel.empty()
    ch_annotation = (params.annotation && params.annotation instanceof String)
        ? channel.fromSamplesheet("annotation", parameters_schema: "${projectDir}/nextflow_schema.json")
        : channel.empty()
    ch_mirna = params.mature && params.mirna_expression ? channel.value([[id: "mirna"], file(params.mirna_expression, checkIfExists: true)]) : channel.empty()

    CIRCRNA(
        ch_samplesheet,
        ch_phenotype,
        ch_fasta,
        ch_gtf,
        ch_blacklist,
        ch_mature,
        ch_annotation,
        ch_versions,
        ch_mirna,
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

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CIRCRNA (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_CIRCRNA.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
