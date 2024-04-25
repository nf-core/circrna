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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CIRCRNA  } from './workflows/circrna'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_circrna_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_circrna_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_circrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta      = getGenomeAttribute('fasta')
params.gtf        = getGenomeAttribute('gtf')
params.bwa        = getGenomeAttribute('bwa')
params.star       = getGenomeAttribute('star')
params.bowtie     = getGenomeAttribute('bowtie')
params.bowtie2    = getGenomeAttribute('bowtie2')
params.mature     = getGenomeAttribute('mature')
params.species_id = getGenomeAttribute('species_id')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CIRCRNA {
    main:

    ch_versions = Channel.empty()

    //
    // WORKFLOW: Run nf-core/circrna workflow
    //
    ch_samplesheet = Channel.fromSamplesheet("input")
    ch_fasta       = Channel.value([[id: "fasta"], file(params.fasta, checkIfExists:true)])
    ch_gtf         = Channel.value([[id: "gtf"], file(params.gtf, checkIfExists:true)])
    ch_mature      = params.module.split(',').contains('mirna_prediction') ? Channel.value([[id: "mature"], file(params.mature, checkIfExists:true)]) : Channel.empty()
    ch_phenotype   = Channel.value([[id: "phenotype"], file(params.phenotype, checkIfExists:true)])
    ch_species     = params.module.split(',').contains('differential_expression') ? Channel.value(params.species_id) : Channel.empty()

    CIRCRNA (
        ch_samplesheet,
        ch_phenotype,
        ch_fasta,
        ch_gtf,
        ch_mature,
        ch_species,
        ch_versions
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
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CIRCRNA ()

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
