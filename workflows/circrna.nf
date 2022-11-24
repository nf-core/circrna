/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCircrna.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Genome params
params.fasta   = params.genome  ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf     = params.genome  ? params.genomes[ params.genome ].gtf ?: false : false
params.bwa     = params.genome && params.tool.contains('ciriquant') ? params.genomes[ params.genome ].bwa ?: false : false
params.star    = params.genome && ( params.tool.contains('circexplorer2') || params.tool.contains('dcc') || params.tool.contains('circrna_finder') ) ? params.genomes[ params.genome ].star ?: false : false
params.bowtie  = params.genome && params.tool.contains('mapsplice') ? params.genomes[ params.genome ].bowtie ?: false : false
params.bowtie2 = params.genome && params.tool.contains('find_circ') ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.mature  = params.genome && params.module.contains('mirna_prediction') ? params.genomes[ params.genome ].mature ?: false : false
params.species = params.genome  ? params.genomes[ params.genome ].species_id ?: false : false

ch_phenotype   = params.phenotype && params.module.contains('differential_expression') ? file(params.phenotype, checkIfExists:true) : Channel.empty()
ch_fasta       = params.fasta ? file(params.fasta) : 'null'
ch_gtf         = params.gtf ? file(params.gtf) : 'null'
ch_mature      = params.mature && params.module.contains('mirna_prediction') ? file(params.mature) : Channel.empty()
ch_species     = params.genome ? Channel.value(params.species) : Channel.value(params.species)


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

// MODULES:

// SUBWORKFLOWS:
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { CIRCRNA_DISCOVERY } from '../subworkflows/local/circrna_discovery'
include { MIRNA_PREDICTION  } from '../subworkflows/local/mirna_prediction'
include { DIFFERENTIAL_EXPRESSION } from '../subworkflows/local/differential_expression'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES:
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'

// SUBWORKFLOWS:
include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CIRCRNA {

    ch_versions = Channel.empty()

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    // Validate input samplesheet & phenotype file
    INPUT_CHECK (
        ch_input,
        ch_phenotype
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .dump(tag: 'map')
    .groupTuple(by: [0])
    .dump(tag: 'group')
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

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
    bowtie2_index  = params.fasta ? params.bowtie2 ? Channel.fromPath(params.bowtie2) : PREPARE_GENOME.out.bowtie2 : []
    bwa_index      = params.fasta ? params.bwa ? Channel.fromPath(params.bwa) : PREPARE_GENOME.out.bwa : []
    chromosomes    = ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) ? PREPARE_GENOME.out.chromosomes : []
    hisat2_index   = params.tool.contains('ciriquant') ? PREPARE_GENOME.out.hisat2 : []
    star_index     = params.fasta ? params.star ? Channel.fromPath(params.star) : PREPARE_GENOME.out.star : []
    segemehl_index = params.fasta ? params.segemehl ? Channel.fromPath(params.segemehl) : PREPARE_GENOME.out.segemehl : []
    ch_versions    = ch_versions.mix(PREPARE_GENOME.out.versions)

    // MODULE: Run FastQC, trimgalore!
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    reads_for_circrna = FASTQC_TRIMGALORE.out.reads
    reads_for_diff_exp = FASTQC_TRIMGALORE.out.reads

    //
    // 2. circRNA Discovery
    //

    CIRCRNA_DISCOVERY(
        reads_for_circrna,
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
        params.duplicates_fun
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

    DIFFERENTIAL_EXPRESSION(
        reads_for_diff_exp,
        ch_gtf,
        ch_fasta,
        hisat2_index,
        PREPARE_GENOME.out.splice_sites
    )

    ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MODULE: MultiQC
    workflow_summary    = WorkflowCircrna.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCircrna.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
