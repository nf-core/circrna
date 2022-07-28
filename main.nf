#!/usr/bin/env nextflow

/*
================================================================================
                            nf-core/circrna
================================================================================
Started August 2020.
Dev version to nf-core Feb 2021.
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/nf-core/circrna
 -------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/circrna
--------------------------------------------------------------------------------
 @Authors
 Barry Digby (@BarryDigby)
--------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

/*
================================================================================
                            Print Help
================================================================================
*/

def json_schema = "$workflow.projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/circrna -profile singularity --input '*_R{1,2}.fastq.gz' --input_type 'fastq' --genome 'GRCh38' --module 'circrna_discovery, mirna_prediction, differential_expression' --tool 'CIRCexplorer2' --phenotype 'metadata.csv' "
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

/*
================================================================================
                        Check parameters
================================================================================
*/

if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)){
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

// Check Tools selected
// Is schema enumerate case sensitive? Seems easier to let toLowerCase take care of misplaced caps from user.
toolList = defineToolList()
tool = params.tool ? params.tool.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tool, toolList)) exit 1, "[nf-core/circrna] error: Unknown tool selected, see --help for more information."

// Check Modules selected
moduleList = defineModuleList()
module = params.module ? params.module.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(module, moduleList)) exit 1, "[nf-core/circrna] error: Unknown module selected, see --help for more information."

// Check phenotype file and stage the channel
// (Must not have NA's, must have 'condition' as colname)

if(params.phenotype){
    pheno_file = file(params.phenotype)
    ch_phenotype = examine_phenotype(pheno_file)
} else {
    ch_phenotype = Channel.empty()
}

// Check BBDUK params
/*
  check adapters file exists
  check combinations of parameters have been supplied correctly
*/

if(params.trim_fastq){
    if(params.adapters){
    adapters = file(params.adapters, checkIfExists: true)
    if(!params.k && !params.ktrim || !params.k && params.ktrim || params.k && !params.ktrim){
        exit 1, "[nf-core/circrna] error: Adapter file provided for trimming but missing values for '--k' and/or '--ktrim'.Please provide values for '--k' and '--ktrim'.\n\nPlease check the parameter documentation online."
    }
    }
    if(params.trimq && !params.qtrim || !params.trimq && params.qtrim){
        exit 1, "[nf-core/circrna] error: Both '--trimq' and '--qtrim' are required to perform quality filtering - only one has been provided.\n\nPlease check the parameter documentation online."
    }
}

// Check filtering params
tools_selected = tool.size()

// Check '--tool_filter'' does not exceed number of tools selected
if(tools_selected > 1 && params.tool_filter > tools_selected){
  exit 1, "[nf-core/circrna] error: The parameter '--tool_filter' (${params.tool_filter}) exceeds the number of tools selected (${params.tool}). Please select a value less than or equal to the number of quantification tools selected ($tools_selected).\n\nPlease check the help documentation."
}

// Check Input data (if !csv choose path)

if(has_extension(params.input, "csv")){

    csv_file = file(params.input, checkIfExists: true)
    ch_input = extract_data(csv_file)

}else if(params.input && !has_extension(params.input, "csv")){

    log.info ""
    log.info "Input data log info:"
    log.info "No input sample CSV file provided, attempting to read from path instead."
    log.info "Reading input data from path: ${params.input}\n"
    log.info ""

    ch_input = retrieve_input_paths(params.input, params.input_type)

}

(bam_input, fastq_input) = ch_input.into(2)

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_multiqc_config = file("$workflow.projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$workflow.projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$workflow.projectDir/docs/images/", checkIfExists: true)

/*
================================================================================
                        PRINTING PARAMETER SUMMARY
================================================================================
*/

log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName

def summary = [:]
if (workflow.revision)        summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = custom_runName ?: workflow.runName
if (workflow.containerEngine) summary['Container'] = "${workflow.containerEngine} - ${workflow.container}"
summary['Max Resources']      = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
summary['Config Files']       = workflow.configFiles.join(', ')
summary['Launch dir']         = workflow.launchDir
summary['Output dir']         = params.outdir
summary['Publish dir mode']   = params.publish_dir_mode
summary['Working dir']        = workflow.workDir
summary['Script dir']         = workflow.projectDir
summary['User']               = workflow.userName

summary['Input']              = params.input
summary['Input type']         = params.input_type
if('differential_expression' in module) summary['Phenotype file'] = params.phenotype
summary['circRNA tool(s)']    = params.tool
summary['modules']            = params.module
summary['BSJ filter']         = params.bsj_reads
if(tools_selected > 1) summary['Tool filter'] = params.tool_filter

summary['Genome']             = params.genome
if(params.fasta)              summary['Reference FASTA']   = params.fasta
if(params.gtf)                summary['Reference GTF']     = params.gtf
if(params.bowtie)             summary['Bowtie indices']    = params.bowtie
if(params.bowtie2)            summary['Bowtie2 indices']   = params.bowtie2
if(params.bwa)                summary['BWA indices']       = params.bwa
if(params.fasta_fai)          summary['SAMtools index']    = params.fasta_fai
if(params.hisat)              summary['HISAT2 indices']    = params.hisat2
if(params.star)               summary ['STAR indices']     = params.star
if(params.segemehl)           summary ['Segemehl index']   = params.segemehl

summary['Save reference files']              = params.save_reference
summary['Save QC intermediates']             = params.save_qc_intermediates
summary['Save RNASeq intermediates']         = params.save_rnaseq_intermediates
summary['Save Quantification intermediates'] = params.save_quantification_intermediates
summary['Save miRNA predictions']            = params.save_mirna_predictions

summary['Skip BBDUK']     = params.trim_fastq
if(params.trim_fastq){
                            summary['BBDUK']             = "Enabled"
if(params.adapters)        summary['Adapter file']      = params.adapters
if(params.k)               summary['k']                 = params.k
if(params.ktrim)           summary['ktrim']             = params.ktrim
if(params.hdist)           summary['hdist']             = params.hdist
if(params.trimq)           summary['trimq']             = params.trimq
if(params.qtrim)           summary['qtrim']             = params.qtrim
if(params.minlen)          summary['minlen']            = params.minlen
}

if('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool){
if(params.alignIntronMax)                      summary['alignIntronMax']               = params.alignIntronMax
if(params.alignIntronMin)                      summary['alignIntronMin']               = params.alignIntronMin
if(params.alignMatesGapMax)                    summary['alignMatesGapMax']             = params.alignMatesGapMax
if(params.alignSJDBoverhangMin)                summary['alignSJDBoverhangMin']         = params.alignSJDBoverhangMin
if(params.alignSJoverhangMin)                  summary['alignSJoverhangMin']           = params.alignSJoverhangMin
if(params.alignSoftClipAtReferenceEnds)        summary['alignSoftClipAtReferenceEnds'] = params.alignSoftClipAtReferenceEnds
if(params.alignTranscriptsPerReadNmax)         summary['alignTranscriptsPerReadNmax']  = params.alignTranscriptsPerReadNmax
if(params.chimJunctionOverhangMin)             summary['chimJunctionOverhangMin']      = params.chimJunctionOverhangMin
if(params.chimScoreMin)                        summary['chimScoreMin']                 = params.chimScoreMin
if(params.chimScoreSeparation)                 summary['chimScoreSeparation']          = params.chimScoreSeparation
if(params.chimSegmentMin)                      summary['chimSegmentMin']               = params.chimSegmentMin
if(params.genomeLoad)                          summary['genomeLoad']                   = params.genomeLoad
if(params.limitSjdbInsertNsj)                  summary['limitSjdbInsertNsj']           = params.limitSjdbInsertNsj
if(params.outFilterMatchNminOverLread)         summary['outFilterMatchNminOverLread']  = params.outFilterMatchNminOverLread
if(params.outFilterMismatchNoverLmax)          summary['outFilterMismatchNoverLmax']   = params.outFilterMismatchNoverLmax
if(params.outFilterMultimapNmax)               summary['outFilterMultimapNmax']        = params.outFilterMultimapNmax
if(params.outFilterMultimapScoreRange)         summary['outFilterMultimapScoreRange']  = params.outFilterMultimapScoreRange
if(params.outFilterScoreMinOverLread)          summary['outFilterScoreMinOverLread']   = params.outFilterScoreMinOverLread
if(params.outSJfilterOverhangMin)              summary['outSJfilterOverhangMin']       = params.outSJfilterOverhangMin
if(params.sjdbOverhang)                        summary['sjdbOverhang']                 = params.sjdbOverhang
if(params.sjdbScore)                           summary['sjdbScore']                    = params.sjdbScore
if(params.winAnchorMultimapNmax)               summary['winAnchorMultimapNmax']        = params.winAnchorMultimapNmax
}

if(workflow.profile.contains('awsbatch')){
    summary['AWS Region'] = params.awsregion
    summary['AWS Queue']  = params.awsqueue
    summary['AWS CLI']    = params.awscli
}

if(params.email || params.email_on_fail){
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-circrna-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/circrna Workflow Summary'
    section_href: 'https://github.com/nf-core/circrna'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


/*
================================================================================
                            Stage Parameters
================================================================================
*/

/* Note to reviewer:
 * I get warning messages when the parameters are not used by the workflow (for e.g when tool that uses bowtie not used):
 * 'WARN: Access to undefined parameter `bowtie` -- Initialise it to a default value eg. `params.bowtie = some_value`'
 * Is there a way to stop these warnings from printing or is it ok to leave as is
 */

params.fasta     = params.genome ? params.genomes[params.genome].fasta ?: false : false
//params.fasta_fai = params.genome ? params.genomes[params.genome].fasta_fai ?: false : false # I don't trust iGenomes has this for all species. Trivial to make in terms of run time.
params.gtf       = params.genome ? params.genomes[params.genome].gtf ?: false : false
params.bwa       = params.genome && 'ciriquant' in tool ? params.genomes[params.genome].bwa ?: false : false
params.star      = params.genome && ('circexplorer2' || 'dcc' || 'circrna_finder' in tool) ? params.genomes[params.genome].star ?: false : false
params.bowtie    = params.genome && 'mapsplice' in tool ? params.genomes[params.genome].bowtie ?: false : false
params.bowtie2   = params.genome && 'find_circ' in tool ? params.genomes[params.genome].bowtie2 ?: false : false
params.mature    = params.genome && 'mirna_prediction' in module ? params.genomes[params.genome].mature ?: false : false
params.species   = params.genome ? params.genomes[params.genome].species_id?: false : false

ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : 'null'
ch_gtf = params.gtf ? Channel.value(file(params.gtf)) : 'null'
ch_mature = params.mature && 'mirna_prediction' in module ? Channel.value(file(params.mature)) : 'null'
ch_species = params.genome ? Channel.value(params.species) : Channel.value(params.species)

/*
================================================================================
                            SOFTWARE VERSIONS
================================================================================
*/

/*
  Note to reviewer:
  I struggled capturing tool versions, particularly with the regexes in 'scrape_software_versions.py'
  For now, the process is hardcoded.
*/

process SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
    saveAs: {filename ->
        if (filename.indexOf(".tsv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.tsv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    echo "37.62" > v_bbduk.txt
    echo "2.29.2" > v_bedtools.txt
    echo "1.2.3" > v_bowtie.txt
    echo "2.3.5.1" > v_bowtie2.txt
    echo "0.7.17" > v_bwa.txt
    echo "2.3.8" > v_circexplorer2.txt
    echo "1.1.1" > v_ciriquant.txt
    echo "8.0.92" > v_java.txt
    echo "2.2.1" > v_mapsplice.txt
    echo "3.3a" > v_miranda.txt
    echo "5.26.2" > v_perl.txt
    echo "3.6.3" > v_R.txt
    echo "0.3.4" > v_segemehl.txt
    echo "2.24.1" > v_picard.txt
    echo "2.7.15" > v_python.txt
    echo "1.10" > v_samtools.txt
    echo "2.6.1d" > v_star.txt
    echo "2.1.1" > v_stringtie.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
================================================================================
                            BUILD INDICES
================================================================================
*/

process BWA_INDEX {
    tag "${fasta}"
    label 'proces_medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/${it}" : null }

    when:
    !params.bwa && !params.genome && params.fasta && 'ciriquant' in tool && 'circrna_discovery' in module

    input:
    file(fasta) from ch_fasta

    output:
    file("BWAIndex") into bwa_built

    script:
    """
    mkdir -p BWAIndex
    bwa index $fasta -p BWAIndex/${fasta.baseName}
    """
}

/*
 * Note to reviewer:
 * I'll verbalise what I think I am doing here, please let me know if there are contradictions in the code.
 * 1. If igenomes '--genome' param passed to script, use pre-built indices from igenomes.
 * 2. If path to indices is provided && --genome is null, use the path.
 * 3. If neither are available, last resort is to build the indices.
 *
 * Had to include '&& ciriquant' below or else it attempts to stage 'false' in file channel, breaks workflow.
 * This does not happen with other index channels which is confusing.
*/

ch_bwa = params.genome && 'ciriquant' in tool ? Channel.value(file(params.bwa)) : params.bwa && 'ciriquant' in tool ? Channel.value(file(params.bwa)) : bwa_built

process SAMTOOLS_INDEX {
    tag "${fasta}"
    label 'process_low'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/SamtoolsIndex/${it}" : null }

    when:
    !params.fasta_fai && params.fasta

    input:
    file(fasta) from ch_fasta

    output:
    file("${fasta}.fai") into fai_built

    script:
    """
    samtools faidx ${fasta}
    """
}

ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fai_built

process HISAT2_INDEX {
    tag "${fasta}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/Hisat2Index/${it}" : null }

    when:
    !params.hisat && params.fasta && ('differential_expression' in module || 'ciriquant' in tool)

    input:
    file(fasta) from ch_fasta

    output:
    file("${fasta.baseName}.*.ht2") into hisat_built
    val("${workflow.launchDir}/${params.outdir}/reference_genome/Hisat2Index") into hisat_path

    script:
    """
    hisat2-build \\
        -p ${task.cpus} \\
        $fasta \\
        ${fasta.baseName}
    """
}

// Hisat2 not available from igenomes so this is more straight forward, path or else build.
// CIRIquant only wants path to index files, does not need files symlinked in work dir hence path not built here.
ch_hisat = params.hisat ? Channel.value(file(params.hisat)) : hisat_path

process STAR_INDEX {
    tag "${fasta}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/${it}" : null }

    when:
    !params.star && params.fasta && params.gtf && ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    input:
    file(fasta) from ch_fasta
    file(gtf) from ch_gtf

    output:
    file("STARIndex") into star_built

    script:
    """
    mkdir -p STARIndex

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbOverhang ${params.sjdbOverhang} \\
        --sjdbGTFfile $gtf \\
        --genomeDir STARIndex/ \\
        --genomeFastaFiles $fasta
    """
}

ch_star = params.genome ? Channel.value(file(params.star)) : params.star ? Channel.value(file(params.star)) : star_built

process BOWTIE_INDEX {
    tag "${fasta}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/BowtieIndex/${it}" : null }

    when:
    !params.bowtie && params.fasta && 'mapsplice' in tool && 'circrna_discovery' in module

    input:
    file(fasta) from ch_fasta

    output:
    file ("${fasta.baseName}.*") into bowtie_built

    script:
    """
    bowtie-build \\
        --threads ${task.cpus} \\
        $fasta \\
        ${fasta.baseName}
    """
}

ch_bowtie = params.genome ? Channel.fromPath("${params.bowtie}*") : params.bowtie ? Channel.fromPath("${params.bowtie}*", checkIfExists: true).ifEmpty { exit 1, "[nf-core/circrna] error: Bowtie index directory not found: ${params.bowtie}"} : bowtie_built

process BOWTIE2_INDEX {
    tag "${fasta}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/Bowtie2Index/${it}" : null }

    when:
    !params.bowtie2 && params.fasta && 'find_circ' in tool && 'circrna_discovery' in module

    input:
    file(fasta) from ch_fasta

    output:
    file ("${fasta.baseName}.*") into bowtie2_built

    script:
    """
    bowtie2-build \\
        --threads ${task.cpus} \\
        $fasta \\
        ${fasta.baseName}
    """
}

ch_bowtie2 = params.genome ? Channel.fromPath("${params.bowtie2}*") : params.bowtie2 ? Channel.fromPath("${params.bowtie2}*", checkIfExists: true).ifEmpty { exit 1, "[nf-core/circrna] error: Bowtie2 index directory not found: ${params.bowtie2}"} : bowtie2_built
(ch_bowtie2_anchors, ch_bowtie2_find_circ) = ch_bowtie2.into(2)

process SEGEMEHL_INDEX{
    tag "${fasta}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/SegemehlIndex/${it}" : null }

    when:
    !params.segemehl && params.fasta && 'segemehl' in tool && 'circrna_discovery' in module

    input:
    file(fasta) from ch_fasta

    output:
    file("${fasta.baseName}.idx") into segemehl_built

    script:
    """
    segemehl.x \\
        -t ${task.cpus} \\
        -d $fasta \\
        -x "${fasta.baseName}.idx"
    """
}

ch_segemehl = params.segemehl ? Channel.fromPath("${params.segemehl}*.idx", checkIfExists: true) : segemehl_built

/*
================================================================================
                            Misc circRNA Requirements
================================================================================
*/

process FILTER_GTF{
    tag"${gtf}"

    when:
    'circrna_discovery' in module

    input:
    file(gtf) from ch_gtf

    output:
    file("filt.gtf") into ch_gtf_filtered

    script:
    """
    grep -vf ${workflow.projectDir}/bin/unwanted_biotypes.txt $gtf > filt.gtf
    """
}


if(('mapsplice' in tool || 'find_circ' in tool) && 'circrna_discovery' in module){
    file("${params.outdir}/reference_genome/chromosomes").mkdirs()
    ch_fasta.splitFasta(record: [id:true])
            .map{ record -> record.id.toString() }
            .flatten()
            .set{ ID }

    ch_fasta.splitFasta(file: true)
            .flatten()
            .merge(ID).map{ it ->
                            file = it[0]
                            chr_id = it[1]
                            file.copyTo("${params.outdir}/reference_genome/chromosomes/${chr_id}.fa")
                          }

    stage_chromosomes = Channel.value("${workflow.launchDir}/${params.outdir}/reference_genome/chromosomes")
}

ch_chromosomes = ('mapsplice' in tool || 'find_circ' in tool) ? stage_chromosomes : 'null'

/*
 * DEBUG
 * No signature of method: nextflow.util.BlankSeparatedList.toRealPath() is applicable for argument types: () values: []
 * error in below YML process (--genome, no params passed to gtf/fasta/etc.. )
 * I suspect this is because bwa_built is a directory (BWAIndex), whilst params.bwa from iGenomes is a collection of files
 * (structure of bwa on iGenomes is frustrating)
 * If iGenomes bwa used, place the files in a directory so it matches bwa_built and "path" method.
 */

if(params.genome && 'ciriquant' in tool){
    file("${params.outdir}/reference_genome/BWAIndex").mkdirs()
    ch_bwa.flatten().map{ it -> it.copyTo("${params.outdir}/reference_genome/BWAIndex")}
    ch_bwa = Channel.value(file("${params.outdir}/reference_genome/BWAIndex"))
}

process CIRIQUANT_YML{

    when:
    'ciriquant' in tool && 'circrna_discovery' in module

    input:
    file(gtf) from ch_gtf
    file(fasta) from ch_fasta
    file(bwa) from ch_bwa
    val(hisat) from ch_hisat

    output:
    file("travis.yml") into ch_ciriquant_yml

    script:
    bwa_prefix = fasta.toString() == 'genome.fa' ? fasta.toString() : fasta.toString() - ~/.(fa|fasta)$/
    hisat_prefix = fasta.toString() - ~/.(fa|fasta)$/
    fasta_path = fasta.toRealPath()
    gtf_path = gtf.toRealPath()
    bwa_path = bwa.toRealPath()

    """
    BWA=`whereis bwa | cut -f2 -d':'`
    HISAT2=`whereis hisat2 | cut -f2 -d':'`
    STRINGTIE=`whereis stringtie | cut -f2 -d':'`
    SAMTOOLS=`whereis samtools | cut -f2 -d':' | awk '{print \$1}'`

    touch travis.yml
    printf "name: ciriquant\n\
    tools:\n\
     bwa: \$BWA\n\
     hisat2: \$HISAT2\n\
     stringtie: \$STRINGTIE\n\
     samtools: \$SAMTOOLS\n\n\
    reference:\n\
     fasta: ${fasta_path}\n\
     gtf: ${gtf_path}\n\
     bwa_index: ${bwa_path}/${bwa_prefix}\n\
     hisat_index: ${hisat}/${hisat_prefix}" >> travis.yml
    """
}

process GENE_ANNOTATION{
    tag "${gtf}"
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/${it}" : null }

    when:
    !params.circexplorer2_annotation && params.gtf && ('circexplorer2' || 'mapsplice' in tool) && 'circrna_discovery' in module

    input:
    file(gtf) from ch_gtf

    output:
    file("${gtf.baseName}.txt") into ch_gene_txt

    script:
    """
    gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} ${gtf.baseName}.genepred
    perl -alne '\$"="\t";print "@F[11,0..9]"' ${gtf.baseName}.genepred > ${gtf.baseName}.txt
    """
}

ch_gene = params.circexplorer2_annotation ? Channel.value(file(params.circexplorer2_annotation)) : ch_gene_txt

/*
================================================================================
                            Stage Input Data
================================================================================
*/

process BAM_TO_FASTQ{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_qc_intermediates ? "quality_control/SamToFastq/${it}" : null }

    when:
    params.input_type == 'bam'

    input:
    tuple val(base), file(bam) from bam_input

    output:
    tuple val(base), file('*.fq.gz') into fastq_built

    script:
    """
    picard \\
        -Xmx${task.memory.toGiga()}g \\
        SamToFastq \\
        I=$bam \\
        F=${base}_R1.fq.gz \\
        F2=${base}_R2.fq.gz \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

// fastq_input staged on line 116
if(params.input_type == 'bam'){
    (fastqc_reads, trimming_reads, raw_reads) = fastq_built.into(3)
}else if(params.input_type == 'fastq'){
    (fastqc_reads, trimming_reads, raw_reads) = fastq_input.into(3)
}

process FASTQC_RAW {
    tag "${base}"
    label 'process_low'
    label 'py3'

    input:
    tuple val(base), file(fastq) from fastqc_reads

    output:
    file("*.{html,zip}") into fastqc_raw

    script:
    """
    fastqc -q $fastq --threads ${task.cpus}
    """
}

/*
================================================================================
                                    BBDUK
================================================================================
*/

process BBDUK {
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.fq.gz",
        saveAs: { params.save_qc_intermediates ? "quality_control/BBDUK/${it}" : null }

    when:
    params.trim_fastq

    input:
    tuple val(base), file(fastq) from trimming_reads
    path adapters from params.adapters

    output:
    tuple val(base), file('*.trim.fq.gz') into trim_reads_ch, fastqc_trim_reads
    file("*BBDUK.txt") into bbduk_stats_ch

    script:
    def adapter = params.adapters ? "ref=${params.adapters}" : ''
    def k = params.k ? "k=${params.k}" : ''
    def ktrim = params.ktrim ? "ktrim=${params.ktrim}" : ''
    def hdist = params.hdist ? "hdist=${params.hdist}" : ''
    def trimq = params.trimq ? "trimq=${params.trimq}" : ''
    def qtrim = params.qtrim ? "qtrim=${params.qtrim}" : ''
    def minlen = params.minlen ? "minlen=${params.minlen}" : ''
    """
    bbduk.sh \\
        -Xmx${task.memory.toGiga()}g \\
        threads=${task.cpus} \\
        in1=${fastq[0]} \\
        in2=${fastq[1]} \\
        out1=${base}_R1.trim.fq.gz \\
        out2=${base}_R2.trim.fq.gz \\
        $adapter \\
        $k \\
        $ktrim \\
        $trimq \\
        $qtrim \\
        $minlen \\
        stats=${base}_BBDUK.txt
    """
}

aligner_reads = params.trim_fastq ? trim_reads_ch : raw_reads

process FASTQC_BBDUK {
    tag "${base}"
    label 'process_low'
    label 'py3'

    when:
    params.trim_fastq

    input:
    tuple val(base), file(fastq) from fastqc_trim_reads

    output:
    file ("*.{html,zip}") into fastqc_trimmed

    script:
    """
    fastqc -q $fastq --threads ${task.cpus}
    """
}

(star_pass1_reads, star_pass2_reads, find_circ_reads, ciriquant_reads, mapsplice_reads, segemehl_reads, dcc_mate1_reads, dcc_mate2_reads, hisat_reads) = aligner_reads.into(9)

/*
================================================================================
                    circRNA quantification + annotation
================================================================================
*/

process CIRIQUANT{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/CIRIquant/intermediates/${it}" : null }

    when:
    'ciriquant' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(fastq) from ciriquant_reads
    file(ciriquant_yml) from ch_ciriquant_yml

    output:
    tuple val(base), file("${base}") into ciriquant_intermediates
    tuple val(base), val("CIRIquant"), file("${base}_ciriquant_circs.bed") into ciriquant_annotated
    tuple val(base), file("${base}_ciriquant.bed") into ciriquant_results

    script:
    """
    CIRIquant \\
        -t ${task.cpus} \\
        -1 ${fastq[0]} \\
        -2 ${fastq[1]} \\
        --config $ciriquant_yml \\
        --no-gene \\
        -o ${base} \\
        -p ${base}

    ## Apply Filtering
    cp ${base}/${base}.gtf .

    ## extract counts (convert float/double to int [no loss of information])
    grep -v "#" ${base}.gtf | awk '{print \$14}' | cut -d '.' -f1 > counts
    grep -v "#" ${base}.gtf | awk -v OFS="\t" '{print \$1,\$4,\$5,\$7}' > ${base}.tmp
    paste ${base}.tmp counts > ${base}_unfilt.bed

    ## filter bsj_reads
    awk '{if(\$5 >= ${params.bsj_reads}) print \$0}' ${base}_unfilt.bed > ${base}_filt.bed
    grep -v '^\$' ${base}_filt.bed > ${base}_ciriquant

    ## correct offset bp position
    awk -v OFS="\t" '{\$2-=1;print}' ${base}_ciriquant > ${base}_ciriquant.bed

    rm ${base}.gtf

    ## Re-work for Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_ciriquant.bed > ${base}_ciriquant_circs.bed
    """
}

process STAR_1PASS{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/STAR/1st_Pass/${it}" : null }

    when:
    ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    input:
    tuple val(base), file(reads) from star_pass1_reads
    file(star_idx) from ch_star

    output:
    file("${base}/*SJ.out.tab") into sjdb_ch
    file("${base}") into star_1st_pass_output

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    mkdir -p ${base}

    STAR \\
        --alignIntronMax ${params.alignIntronMax} \\
        --alignIntronMin ${params.alignIntronMin} \\
        --alignMatesGapMax ${params.alignMatesGapMax} \\
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \\
        --alignSJoverhangMin ${params.alignSJoverhangMin} \\
        --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \\
        --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \\
        --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \\
        --chimOutType Junctions SeparateSAMold \\
        --chimScoreMin ${params.chimScoreMin} \\
        --chimScoreSeparation ${params.chimScoreSeparation} \\
        --chimSegmentMin ${params.chimSegmentMin} \\
        --genomeDir ${star_idx} \\
        --genomeLoad ${params.genomeLoad} \\
        --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \\
        --outFileNamePrefix ${base}/${base}. \\
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \\
        --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \\
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \\
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \\
        --outFilterType BySJout \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \\
        ${readFilesCommand} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --sjdbScore ${params.sjdbScore} \\
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

process SJDB_FILE{
    tag "${sjdb}"
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/STAR/SJFile/${it}" : null }

    when:
    ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    input:
    file(sjdb) from sjdb_ch

    output:
    file("*SJFile.tab") into sjdbfile_ch

    shell:
    '''
    base=$(basename !{sjdb} .SJ.out.tab)
    awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' !{sjdb} > ${base}.SJFile.tab
    '''
}

(sjdbfile_pass2, sjdbfile_mate1, sjdbfile_mate2) = sjdbfile_ch.into(3)

process STAR_2PASS{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/STAR/2nd_Pass/${it}" : null }

    when:
    ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    input:
    tuple val(base), file(reads) from star_pass2_reads
    file(sjdbfile) from sjdbfile_pass2.collect()
    file(star_idx) from ch_star

    output:
    tuple val(base), file("${base}/${base}.Chimeric.out.junction") into circexplorer2_input
    tuple val(base), file("${base}") into circrna_finder_input, dcc_pairs


    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    mkdir -p ${base}

    STAR \\
        --alignIntronMax ${params.alignIntronMax} \\
        --alignIntronMin ${params.alignIntronMin} \\
        --alignMatesGapMax ${params.alignMatesGapMax} \\
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \\
        --alignSJoverhangMin ${params.alignSJoverhangMin} \\
        --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \\
        --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \\
        --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \\
        --chimOutType Junctions SeparateSAMold \\
        --chimScoreMin ${params.chimScoreMin} \\
        --chimScoreSeparation ${params.chimScoreSeparation} \\
        --chimSegmentMin ${params.chimSegmentMin} \\
        --genomeDir ${star_idx} \\
        --genomeLoad ${params.genomeLoad} \\
        --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \\
        --outFileNamePrefix ${base}/${base}. \\
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \\
        --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \\
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \\
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \\
        --outFilterType BySJout \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outBAMsortingBinsN 150 \\
        --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \\
        ${readFilesCommand} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --sjdbFileChrStartEnd ${sjdbfile} \\
        --sjdbScore ${params.sjdbScore} \\
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

process CIRCEXPLORER2{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/CIRCexplorer2/intermediates/${it}" : null }

    when:
    'circexplorer2' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(chimeric_reads) from circexplorer2_input
    file(fasta) from ch_fasta
    file(gene_annotation) from ch_gene

    output:
    tuple val(base), file("${base}") into circexplorer2_intermediates
    tuple val(base), file("${base}_circexplorer2.bed") into circexplorer2_results
    tuple val(base), val("CIRCexplorer2"), file("${base}_circexplorer2_circs.bed") into circexplorer2_annotated

    script:
    """
    mkdir -p ${base}

    CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}/${base}.STAR.junction.bed

    CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}/${base}.STAR.junction.bed -o ${base}/${base}.txt

    awk '{if(\$13 >= ${params.bsj_reads}) print \$0}' ${base}/${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_circexplorer2.bed

    ## Re-work for Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_circexplorer2.bed > ${base}_circexplorer2_circs.bed
    """
}

process CIRCRNA_FINDER{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/circRNA_Finder/intermediates/${it}" : null }

    when:
    'circrna_finder' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(star_dir) from circrna_finder_input
    file(fasta) from ch_fasta

    output:
    tuple val(base), file("${base}_circrna_finder.bed") into circrna_finder_results
    tuple val(base), file("${base}") into circrna_finder_intermediates
    tuple val(base), val("circRNA_Finder"), file("${base}_circrna_finder_circs.bed") into circrna_finder_annotated

    script:
    """
    postProcessStarAlignment.pl --starDir ${star_dir}/ --outDir ./

    awk '{if(\$5 >= ${params.bsj_reads}) print \$0}' ${base}.filteredJunctions.bed | awk  -v OFS="\t" -F"\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_circrna_finder.bed

    mkdir -p ${base}

    mv *filteredJunctions* ${base}
    mv *.Chimeric.out.sorted.* ${base}

    ## Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_circrna_finder.bed > ${base}_circrna_finder_circs.bed
    """
}

process DCC_MATE1{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "mate1",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/DCC/intermediates/${base}/${it}" : null }

    when:
    'dcc' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(reads) from dcc_mate1_reads
    file(sjdbfile) from sjdbfile_mate1.collect()
    file(star_idx) from ch_star

    output:
    tuple val(base), file("mate1") into dcc_mate1

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    mkdir -p mate1

    STAR \\
        --alignIntronMax ${params.alignIntronMax} \\
        --alignIntronMin ${params.alignIntronMin} \\
        --alignMatesGapMax ${params.alignMatesGapMax} \\
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \\
        --alignSJoverhangMin ${params.alignSJoverhangMin} \\
        --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \\
        --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \\
        --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \\
        --chimOutType Junctions SeparateSAMold \\
        --chimScoreMin ${params.chimScoreMin} \\
        --chimScoreSeparation ${params.chimScoreSeparation} \\
        --chimSegmentMin ${params.chimSegmentMin} \\
        --genomeDir ${star_idx} \\
        --genomeLoad ${params.genomeLoad} \\
        --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \\
        --outFileNamePrefix mate1/${base}. \\
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \\
        --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \\
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \\
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \\
        --outFilterType BySJout \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outBAMsortingBinsN 150 \\
        --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \\
        ${readFilesCommand} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --sjdbFileChrStartEnd ${sjdbfile} \\
        --sjdbScore ${params.sjdbScore} \\
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

process DCC_MATE2{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "mate2",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/DCC/intermediates/${base}/${it}" : null }

    when:
    'dcc' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(reads) from dcc_mate2_reads
    file(sjdbfile) from sjdbfile_mate2.collect()
    file(star_idx) from ch_star

    output:
    tuple val(base), file("mate2") into dcc_mate2

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    mkdir -p mate2

    STAR \\
        --alignIntronMax ${params.alignIntronMax} \\
        --alignIntronMin ${params.alignIntronMin} \\
        --alignMatesGapMax ${params.alignMatesGapMax} \\
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \\
        --alignSJoverhangMin ${params.alignSJoverhangMin} \\
        --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \\
        --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \\
        --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \\
        --chimOutType Junctions SeparateSAMold \\
        --chimScoreMin ${params.chimScoreMin} \\
        --chimScoreSeparation ${params.chimScoreSeparation} \\
        --chimSegmentMin ${params.chimSegmentMin} \\
        --genomeDir ${star_idx} \\
        --genomeLoad ${params.genomeLoad} \\
        --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \\
        --outFileNamePrefix mate2/${base}. \\
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \\
        --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \\
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \\
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \\
        --outFilterType BySJout \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outBAMsortingBinsN 150 \\
        --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \\
        ${readFilesCommand} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --sjdbFileChrStartEnd ${sjdbfile} \\
        --sjdbScore ${params.sjdbScore} \\
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

ch_dcc_dirs = dcc_pairs.join(dcc_mate1).join(dcc_mate2)

process DCC{
    tag "${base}"
    label 'py3'
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/DCC/intermediates/${base}/${it}" : null }

    when:
    'dcc' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(pairs), file(mate1), file(mate2) from ch_dcc_dirs
    file(gtf) from ch_gtf
    file(fasta) from ch_fasta

    output:
    tuple val(base), file("${base}_dcc.bed") into dcc_results
    tuple val(base), file("${base}") into dcc_intermediates
    tuple val(base), val("DCC"), file("${base}_dcc_circs.bed") into dcc_annotated

    script:
    COJ="Chimeric.out.junction"
    """
    sed -i 's/^chr//g' $gtf
    printf "${base}/${base}.${COJ}" > samplesheet
    printf "mate1/${base}.${COJ}" > mate1file
    printf "mate2/${base}.${COJ}" > mate2file
    DCC @samplesheet -mt1 @mate1file -mt2 @mate2file -D -an $gtf -Pi -ss -F -M -Nr 1 1 -fg -A $fasta -N -T ${task.cpus}

    ## Add strand to counts
    awk '{print \$6}' CircCoordinates >> strand
    paste CircRNACount strand | tail -n +2 | awk -v OFS="\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${base}_dcc.txt

    ## filter reads
    awk '{if(\$5 >= ${params.bsj_reads}) print \$0}' ${base}_dcc.txt > ${base}_dcc.filtered

    ## fix start position (+1)
    awk -v OFS="\t" '{\$2-=1;print}' ${base}_dcc.filtered > ${base}_dcc.bed

    mkdir -p ${base}
    rm strand
    rm ${base}_dcc.txt
    rm ${base}_dcc.filtered
    find . -maxdepth 1 -mindepth 1 -type f -not -name ${base}_dcc.bed -print0 | xargs -0 mv -t ${base}/

    ## Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_dcc.bed > ${base}_dcc_circs.bed
    """
}

process FIND_ANCHORS{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "{*.*}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/find_circ/intermediates/${base}/${it}" : null }

    when:
    'find_circ' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(fastq) from find_circ_reads
    val(fasta) from ch_fasta
    file(bowtie2_index) from ch_bowtie2_anchors.collect()

    output:
    tuple val(base), file("${base}_anchors.qfa.gz") into ch_anchors
    tuple val(base), file("*.*") into find_anchors_intermediates

    script:
    """
    bowtie2 -p ${task.cpus} --very-sensitive --mm -D 20 --score-min=C,-15,0 \\
    -x ${fasta.baseName} -q -1 ${fastq[0]} -2 ${fastq[1]} \\
    | samtools view -hbuS - | samtools sort --threads ${task.cpus} -m 2G - > ${base}.bam

    samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

    unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz
    """
}

// avoid input collision here because iGenomes index comes pre-packaged with ref fasta file.
// Filter out the ref fasta file from pre-built index so I can use file(fasta) when --genome null
ch_avoid_collisions = ch_bowtie2_find_circ.flatten().filter{ file -> file.getFileName().toString().endsWith(".bt2") }

process FIND_CIRC{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.sites.*",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/find_circ/intermediates/${base}/${it}" : null }

    when:
    'find_circ' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(anchors) from ch_anchors
    file(fasta) from ch_fasta
    file(bowtie2_index) from ch_avoid_collisions.collect()
    val(fasta_chr_path) from ch_chromosomes

    output:
    tuple val(base), file("${base}_find_circ.bed") into find_circ_results
    tuple val(base), file("*.sites.*") into find_circ_intermediates
    tuple val(base), val("find_circ"), file("${base}_find_circ_circs.bed") into find_circ_annotated

    script:
    """
    bowtie2 -p ${task.cpus} --reorder --mm -D 20 --score-min=C,-15,0 -q -x ${fasta.baseName} \\
    -U $anchors | python ${workflow.projectDir}/bin/find_circ.py -G $fasta_chr_path -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

    ## filtering
    grep circ ${base}.sites.bed | grep -v chrM | python ${workflow.projectDir}/bin/sum.py -2,3 | python ${workflow.projectDir}/bin/scorethresh.py -16 1 | python ${workflow.projectDir}/bin/scorethresh.py -15 2 | python ${workflow.projectDir}/bin/scorethresh.py -14 2 | python ${workflow.projectDir}/bin/scorethresh.py 7 ${params.bsj_reads} | python ${workflow.projectDir}/bin/scorethresh.py 8,9 35 | python ${workflow.projectDir}/bin/scorethresh.py -17 100000 >> ${base}.txt

    tail -n +2 ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_find_circ.bed

    ## Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_find_circ.bed > ${base}_find_circ_circs.bed
    """
}

process MAPSPLICE_ALIGN{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/MapSplice/intermediates/${it}" : null }

    when:
    'mapsplice' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(fastq) from mapsplice_reads
    val(mapsplice_ref) from ch_chromosomes
    file(bowtie_index) from ch_bowtie.collect()
    file(gtf) from ch_gtf

    output:
    tuple val(base), file("${base}/fusions_raw.txt") into mapsplice_fusion
    tuple val(base), file("${base}") into mapsplice_align_intermediates

    script:
    def prefix = gtf.toString() - ~/.gtf/
    def handleGzip_R1 = fastq[0].toString().endsWith('.gz') ? "gzip -d --force ${fastq[0]}" : ''
    def handleGzip_R2 = fastq[1].toString().endsWith('.gz') ? "gzip -d --force ${fastq[1]}" : ''
    def read1 = fastq[0].toString().endsWith('.gz') ? fastq[0].toString() - ~/.gz/ : fastq[0]
    def read2 = fastq[1].toString().endsWith('.gz') ? fastq[1].toString() - ~/.gz/ : fastq[1]
    """
    $handleGzip_R1
    $handleGzip_R2

    mapsplice.py \\
        -c $mapsplice_ref \\
        -x $prefix \\
        -1 ${read1} \\
        -2 ${read2} \\
        -p ${task.cpus} \\
        --bam \\
        --seglen 25 \\
        --min-intron ${params.alignIntronMin} \\
        --max-intron ${params.alignIntronMax} \\
        --min-map-len 40 \\
        --fusion-non-canonical \\
        --min-fusion-distance 200 \\
        --gene-gtf $gtf \\
        -o $base
    """
}

/*
  STEP 6.7.2: MapSplice quantification
*/

process MAPSPLICE_PARSE{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}/*",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/MapSplice/intermediates/${it}" : null }

    when:
    'mapsplice' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(raw_fusion) from mapsplice_fusion
    file(fasta) from ch_fasta
    file(gene_annotation) from ch_gene

    output:
    tuple val(base), file("${base}_mapsplice.bed") into mapsplice_results
    tuple val(base), file("${base}/*") into mapsplice_intermediates
    tuple val(base), val("MapSplice"), file("${base}_mapsplice_circs.bed") into mapsplice_annotated

    script:
    """
    mkdir -p ${base}

    CIRCexplorer2 parse -t MapSplice $raw_fusion -b ${base}/${base}.mapsplice.junction.bed

    CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}/${base}.mapsplice.junction.bed -o ${base}/${base}.txt

    awk '{if(\$13 >= ${params.bsj_reads}) print \$0}' ${base}/${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_mapsplice.bed

    ## Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_mapsplice.bed > ${base}_mapsplice_circs.bed
    """
}

process SEGEMEHL_ALIGN{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/Segemehl/intermediates/${it}" : null }

    when:
    'segemehl' in tool && 'circrna_discovery' in module

    input:
    tuple val(base), file(fastq) from segemehl_reads
    file(fasta) from ch_fasta
    file(idx) from ch_segemehl

    output:
    tuple val(base), file("${base}") into segemehl_intermediates
    tuple val(base), file("${base}_segemehl.bed") into segemehl_results
    tuple val(base), val("Segemehl"), file("${base}_segemehl_circs.bed") into segemehl_annotated

    script:
    def handleSam = params.save_quantification_intermediates ? "samtools view -hbS ${base}/${base}.sam > ${base}/${base}.bam && rm ${base}/${base}.sam" : "rm -rf ${base}/${base}.sam"
    """
    mkdir -p ${base}

    segemehl.x \\
        -t ${task.cpus} \\
        -d $fasta \\
        -i $idx \\
        -q ${fastq[0]} \\
        -p ${fastq[1]} \\
        -S \\
        -o ${base}/${base}.sam

    $handleSam

    ## Segemehl does not preserve strand information, nor account for it
    ## when collapsing and counting reads using haarz.x. This is my own fix which does.
    grep ';C;' ${base}/${base}.sngl.bed | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6}' | sort | uniq -c | awk -v OFS="\t" '{print \$2,\$3,\$4,\$5,\$1}' > ${base}/${base}_collapsed.bed

    ## Let user filter by BSJ read count param.
    awk -v OFS="\t" -v BSJ=${params.bsj_reads} '{if(\$5>=BSJ) print \$0}' ${base}/${base}_collapsed.bed > ${base}_segemehl.bed

    ## Re-work for Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_segemehl.bed > ${base}_segemehl_circs.bed
    """
}


// collect circrna info from each alignment process
ch_circs = ciriquant_annotated.mix(circexplorer2_annotated, dcc_annotated, circrna_finder_annotated, find_circ_annotated, mapsplice_annotated, segemehl_annotated)


process ANNOTATION{
    tag "${base}:${tool}"
    label 'process_high'
    publishDir "${params.outdir}/circrna_discovery/${tool}/${base}", mode: params.publish_dir_mode, pattern: "${base}.bed"
    publishDir "${params.outdir}/circrna_discovery/${tool}/annotation_logs", mode: params.publish_dir_mode, pattern: "${base}.log"

    input:
    tuple val(base), val(tool), file(circs) from ch_circs
    file(gtf_filt) from ch_gtf_filtered

    output:
    tuple val(base), val(tool), file("${base}.bed") into ch_annotation
    tuple val(base), file("${base}.log") into annotation_logs

    script:
    """
    mv $circs circs.bed
    bash ${workflow.projectDir}/bin/annotate_outputs.sh &> ${base}.log
    mv master_bed12.bed ${base}.bed.tmp

    ## isolate exon blocks
    awk -FS="\t" '{print \$11}' ${base}.bed.tmp > mature_len.tmp

    ## sum exon block values
    awk -v FS="," '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' mature_len.tmp > mature_length

    ## concat to annotation file.
    paste ${base}.bed.tmp mature_length > ${base}.bed
    """

}


process FASTA{
    tag "${base}:${tool}"
    label 'process_high'
    publishDir "${params.outdir}/circrna_discovery/${tool}/${base}", mode: params.publish_dir_mode, pattern: "fasta/*"

    input:
    tuple val(base), val(tool), file(bed) from ch_annotation
    file(fasta) from ch_fasta

    output:
    tuple val(base), file("fasta/*") into ch_mature_len_fasta

    script:
    """
    ## FASTA sequences (bedtools does not like the extra annotation info - split will not work properly)
    cut -d\$'\t' -f1-12 ${base}.bed > bed12.tmp
    bedtools getfasta -fi $fasta -bed bed12.tmp -s -split -name > circ_seq.tmp
    ## clean fasta header
    grep -A 1 '>' circ_seq.tmp | cut -d: -f1,2,3 > circ_seq.fa && rm circ_seq.tmp
    ## output to dir
    mkdir -p fasta
    awk -F '>' '/^>/ {F=sprintf("fasta/%s.fa",\$2); print > F;next;} {print >> F;}' < circ_seq.fa
    """
}


/*
================================================================================
                        Generate circRNA count matrix
================================================================================
*/

quantification_results = ciriquant_results.mix(circexplorer2_results, circrna_finder_results, dcc_results, find_circ_results, mapsplice_results, segemehl_results)

if(tools_selected > 1){
    process MERGE_TOOLS{
        tag "${base}"

        when:
        'circrna_discovery' in module

        input:
        tuple val(base), file(bed) from quantification_results.groupTuple()

        output:
        file("${base}.bed") into sample_counts

        script:
        """
        ## make list of files for R to read
        ls *.bed > samples.csv

        ## Add catch for empty bed file and delete
        bash ${workflow.projectDir}/bin/check_empty.sh

        ## Use intersection of "n" (params.tool_filter) circRNAs called by tools
        ## remove duplicate IDs, keep highest count.
        Rscript ${workflow.projectDir}/bin/consolidate_algorithms_intersection.R samples.csv $params.tool_filter

        mv combined_counts.bed ${base}.bed
        """
    }
    process COUNT_MATRIX_COMBINED{
        publishDir "${params.outdir}/circrna_discovery", pattern: "count_matrix.txt", mode: params.publish_dir_mode

        when:
        'circrna_discovery' in module

        input:
        file(bed) from sample_counts.collect()

        output:
        file("circRNA_matrix.txt") into circRNA_counts
        file("count_matrix.txt") into matrix

        script:
        """
        python ${workflow.projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
        Rscript ${workflow.projectDir}/bin/reformat_count_matrix.R
        """
    }
}else{
    process COUNT_MATRIX_SINGLE{
        publishDir "${params.outdir}/circrna_discovery", pattern: "count_matrix.txt", mode: params.publish_dir_mode

        when:
        'circrna_discovery' in module

        input:
        file(bed) from quantification_results.collect()
        val(tool) from params.tool

        output:
        file("circRNA_matrix.txt") into circRNA_counts
        file("count_matrix.txt") into matrix

        script:
        """
        # Strip tool name from BED files (no consolidation prior to this step for 1 tool)
        for b in *.bed; do
            basename=\${b%".bed"};
            sample_name=\${basename%"_${tool}"};
            mv \$b \${sample_name}.bed
        done

        python ${workflow.projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
        Rscript ${workflow.projectDir}/bin/reformat_count_matrix.R
        """
    }
}

/*
================================================================================
                            miRNA Prediction
================================================================================
*/

process TARGETSCAN_DATABASE{
    when:
    'mirna_prediction' in module

    input:
    file(mature) from ch_mature

    output:
    file("mature.txt") into ch_mature_txt

    script:
    """
    bash ${workflow.projectDir}/bin/targetscan_format.sh $mature
    """
}

//mirna_input = ciriquant_fasta.mix(circexplorer2_fasta, circrna_finder_fasta, dcc_fasta, mapsplice_fasta, find_circ_fasta, segemehl_fasta).unique().transpose()
mirna_input = ch_mature_len_fasta.unique().transpose()

process MIRNA_PREDICTION{
    tag "${base}"
    label 'process_low'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.miRanda.txt",
        saveAs: { params.save_mirna_predictions ? "mirna_prediction/miRanda/${base}/${it}" : null }
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.targetscan.txt",
        saveAs: { params.save_mirna_predictions ? "mirna_prediction/TargetScan/${base}/${it}" : null }
    when:
    'mirna_prediction' in module

    input:
    tuple val(base), file(fasta) from mirna_input
    file(mirbase) from ch_mature
    file(mirbase_txt) from ch_mature_txt

    output:
    tuple val(base), file("*.miRanda.txt"), file("*.targetscan.txt") into mirna_prediction

    script:
    prefix = fasta.toString() - ~/.fa/
    """
    miranda $mirbase $fasta -strict -out ${prefix}.bindsites.out -quiet
    echo "miRNA Target Score Energy_KcalMol Query_Start Query_End Subject_Start Subject_End Aln_len Subject_Identity Query_Identity" | tr ' ' '\t' > ${prefix}.miRanda.txt

    # Add catch here for non hits (supply NAs to outfile)
    # Making the decision that if miRanda fails, then the miRNA analysis for this circRNA exits cleanly.
    # Happy to rework in the future, but do not want pipeline failing on low confidence circRNA calls.
    ## exit code 1 = fail, 0 = success
    if grep -A 1 -q "Scores for this hit:" ${prefix}.bindsites.out;
    then
        grep -A 1 "Scores for this hit:" ${prefix}.bindsites.out | sort | grep ">" | cut -c 2- | tr ' ' '\t' >> ${prefix}.miRanda.txt

        ##format for targetscan
        cat $fasta | grep ">" | sed 's/>//g' > id
        cat $fasta | grep -v ">" > seq
        echo "0000" > species
        paste id species seq > ${prefix}_ts.txt

        # run targetscan
        targetscan_70.pl mature.txt ${prefix}_ts.txt ${prefix}.targetscan.txt
    else
        ## Add NA's to miRanda cols:
        printf "%0.sNA\t" {1..11} >> ${prefix}.miRanda.txt
        ## Construct TargetScan header
        echo "a_Gene_ID miRNA_family_ID species_ID MSA_start MSA_end UTR_start UTR_end Group_num Site_type miRNA_in_this_species Group_type Species_in_this_group Species_in_this_group_with_this_site_type ORF_overlap" | tr ' ' '\t' > ${prefix}.targetscan.txt
        ## Add NA's to file
        printf "%0.sNA\t" {1..13} >> ${prefix}.targetscan.txt
    fi
    """
}

process MIRNA_TARGETS{
    tag "${base}"
    label 'process_low'
    publishDir "${params.outdir}/mirna_prediction/${base}", mode: params.publish_dir_mode, pattern: "*miRNA_targets.txt"
    publishDir "${params.outdir}/mirna_prediction/${base}/pdf", mode: params.publish_dir_mode, pattern: "*.pdf"

    input:
    tuple val(base), file(miranda), file(targetscan) from mirna_prediction
    file(fasta) from ch_fasta
    file(fai) from ch_fai
    file(filt_gtf) from ch_gtf_filtered
    val(species) from ch_species

    output:
    tuple val(base), file("*.pdf") into circos_plots
    tuple val(base), file("*miRNA_targets.txt") into circrna_mirna_targets

    script:
    def species_id = species + "-"
    """
    ## As before, we have a catch for NA miRNA pred files.
    grep -v "miRNA" $miranda | if grep -q "NA";
    then
        touch ${base}_fail_catch_miRNA_targets.txt
        touch ${base}_fail_catch.pdf
    else
        ## use file name to derive bed12 coordiantes.
        echo *.miRanda.txt | sed -E 's/^(chr[^:]+):([0-9]+)-([0-9]+):([^.]+).*/\\1\\t\\2\\t\\3\\t\\4/' | awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, "0", \$4}' > circs.bed
        bash ${workflow.projectDir}/bin/annotate_outputs.sh &> circ.log
        mv master_bed12.bed circ.bed.tmp

        ## Prep exon track for circlize
        cut -d\$'\t' -f1-12 circ.bed.tmp > bed12.tmp
        bash ${workflow.projectDir}/bin/prep_circos.sh bed12.tmp

        ## add mature spl len (+ 1 bp)
        awk '{print \$11}' circ.bed.tmp | awk -F',' '{for(i=1;i<=NF;i++) printf "%s\\n", \$i}' | awk 'BEGIN {total=0} {total += \$1} END {print total + 1}' > ml
        paste circ.bed.tmp ml > circ.bed

        Rscript ${workflow.projectDir}/bin/mirna_circos.R circ.bed $miranda $targetscan circlize_exons.txt $species_id
    fi
    """
}

/*
================================================================================
                            Differential Expression
================================================================================
*/

// converseley to CIRIquant using Hisat2, this proc does need the index files. stage as new channel
ch_hisat_index_files = params.hisat ? Channel.value(file("${params.hisat}/*")) :  hisat_built

process HISAT_ALIGN{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}.bam",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/Hisat2/${it}" : null }

    when:
    'differential_expression' in module

    input:
    tuple val(base), file(fastq) from hisat_reads
    file(hisat2_index) from ch_hisat_index_files.collect()
    file(fasta) from ch_fasta

    output:
    tuple val(base), file("${base}.bam") into hisat_bam

    script:
    if(fastq[1]){
        """
        hisat2 -p ${task.cpus} --dta -q -x ${fasta.baseName} -1 ${fastq[0]} -2 ${fastq[1]} -t | samtools view -bS - | samtools sort --threads ${task.cpus} -m 2G - > ${base}.bam
        """
    }
    else{
        """
        hisat2 -p ${task.cpus} --dta -q -x ${fasta.baseName} -U ${fastq[0]} -t | samtools view -bS - | samtools sort --threads ${task.cpus} -m 2G - > ${base}.bam
        """
    }
}

process STRINGTIE{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/StringTie/${it}" : null }

    when:
    'differential_expression' in module

    input:
    tuple val(base), file(bam) from hisat_bam
    file(gtf) from ch_gtf

    output:
    file("${base}") into stringtie_dir

    script:
    """
    mkdir ${base}/
    stringtie $bam -e -G $gtf -C ${base}/${base}_cov.gtf -p ${task.cpus} -o ${base}/${base}.gtf -A ${base}/${base}_genes.list
    """
}

process DEA{
    label 'process_medium'
    publishDir "${params.outdir}/differential_expression", pattern: "circRNA", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "boxplots", mode: params.publish_dir_mode
    publishDir "${params.outdir}/quality_control", pattern: "DESeq2_QC", mode: params.publish_dir_mode
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "RNA-Seq",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/${it}" : null }

    when:
    'differential_expression' in module

    input:
    file(gtf_dir) from stringtie_dir.collect()
    file(circ_matrix) from circRNA_counts
    file(phenotype) from ch_phenotype
    val(species) from ch_species

    output:
    file("RNA-Seq") into rnaseq_dir
    file("circRNA") into circrna_dir
    file("boxplots") into boxplots_dir
    file("DESeq2_QC") into qc_plots

    script:
    """
    for i in \$(ls -d */); do sample=\${i%"/"}; file=\${sample}.gtf; touch samples.txt; printf "\$sample\t\${i}\${file}\n" >> samples.txt; done

    prepDE.py -i samples.txt

    ## prepDE && circRNA counts headers are sorted where uppercase preceedes lowercase i.e Z before a
    ## reformat the phenotype file to match the order of the samples.
    head -n 1 $phenotype > header
    tail -n +2 $phenotype | LC_COLLATE=C sort > sorted_pheno
    cat header sorted_pheno > tmp && rm phenotype.csv && mv tmp phenotype.csv

    Rscript ${workflow.projectDir}/bin/DEA.R gene_count_matrix.csv $phenotype $circ_matrix $species ${workflow.projectDir}/bin/ensemblDatabase_map.txt

    mv gene_count_matrix.csv RNA-Seq
    mv transcript_count_matrix.csv RNA-Seq
    """
}

/*
================================================================================
                                    MultiQC
================================================================================
*/

process MULTIQC{
    label 'py3'
    label 'process_low'
    publishDir "${params.outdir}/quality_control/MultiQC", mode: params.publish_dir_mode,
        pattern: "*.html"

    input:
    file(raw_fastqc) from fastqc_raw.collect().ifEmpty([])
    file(trim_fastqc) from fastqc_trimmed.collect().ifEmpty([])
    file(BBDUK_stats) from bbduk_stats_ch.collect().ifEmpty([])
    file(multiqc_config) from ch_multiqc_config
    file(mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file("*.html") into multiqc_out

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file
    """
}

/*
================================================================================
                            Auxiliary functions
================================================================================
*/

// Check integer
def isValidInteger(value){
    value instanceof Integer
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Define list of available tools
def defineToolList() {
    return [
        'ciriquant',
        'circexplorer2',
        'find_circ',
        'circrna_finder',
        'dcc',
        'mapsplice',
        'segemehl'
        ]
}

// Define module list
def defineModuleList() {
    return [
    'circrna_discovery',
    'mirna_prediction',
    'differential_expression'
    ]
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "[nf-core/circrna] error:  Invalid CSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "[nf-core/circrna] error: Cannot find supplied FASTQ or BAM input file. If using input method CSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}"
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Read input files from input CSV
def extract_data(csvFile){
    Channel
        .fromPath(csvFile)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_keys = ["Sample_ID", "Read1", "Read2", "Bam"]
        if(!row.keySet().containsAll(expected_keys)) exit 1, "[nf-core/circrna] error: Invalid CSV input - malformed column names. Please use the column names 'Sample_ID', 'Read1', 'Read2', 'Bam'."

        checkNumberOfItem(row, 4)

        def samples = row.Sample_ID
        def read1 = row.Read1.matches('NA') ? 'NA' : return_file(row.Read1)
        def read2 = row.Read2.matches('NA') ? 'NA' : return_file(row.Read2)
        def bam = row.Bam.matches('NA') ? 'NA' : return_file(row.Bam)

        if(samples == '' || read1 == '' || read2 == '' || bam == '') exit 1, "[nf-core/circrna] error: a field does not contain any information. Please check your CSV file"
        if(read1.matches('NA') && read2.matches('NA') && bam.matches('NA')) exit 1, "[nf-core/circrna] error: A row in your CSV file appears to have missing information."
        if( !read1.matches('NA') && !has_extension(read1, "fastq.gz") && !has_extension(read1, "fq.gz") && !has_extension(read1, "fastq") && !has_extension(read1, "fq")) exit 1, "[nf-core/circrna] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
        if( !read2.matches('NA') && !has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz") && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "[nf-core/circrna] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
        if( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified BAM file either has a non-recognizable extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

        // output tuple mimicking fromFilePairs if FASTQ provided, else tuple for BAM
        if(bam.matches('NA')){
            if(read2.matches('NA')){
                [ samples, read1 ]
            }else{
                [ samples, [read1, read2] ]
            }
        }else{
            [ samples, bam ]
        }

        }
}

// If no input CSV provided, parse input directory containing files.
def retrieve_input_paths(input, type){

    if(type == 'fastq'){

        fastq_files = input
        Channel
                .fromFilePairs(fastq_files)
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty{exit 1, "[nf-core/circrna] error: Your FASTQ files do not have the appropriate extension of either '.fastq.gz', '.fq.gz', .fastq' or '.fq'."}
                .map{ row -> [ row[0], [ row[1][0], row[1][1] ]]}
                .ifEmpty{exit 1, "[nf-core/circrna] error: --input was empty - no files supplied"}
                .set{reads_for_csv}

    }else if(type == 'bam'){

        bam_files = input
        Channel
                .fromFilePairs(bam_files, size: 1)
                .filter{ it =~/.*.bam/}
                .map{ row -> [row[0], [row[1][0]]]}
                .ifEmpty{exit 1, "[nf-core/circrna] error: Cannot find bam file matching: ${bam_files}"}
                .set{reads_for_csv}
    }

    reads_for_csv
                .map{

                def samples = it[0]
                def read1 = (type == 'bam') ? 'NA' : return_file(it[1][0])
                def read2 = (type == 'bam') ? 'NA' : return_file(it[1][1])
                def bam =   (type == 'fastq') ? 'NA' : return_file(it[1][0])

                if(bam.matches('NA')){
                    [ samples, [read1, read2] ]
                }else{
                    [ samples, bam ]
                }

                }
                .ifEmpty{exit 1, "[nf-core/circrna] error: Invalid file paths with --input"}
}


// Check input phenotype file

def examine_phenotype(pheno){

  Channel
        .fromPath(pheno)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_cols = ['condition']

        if (!row.keySet().containsAll(expected_cols)) exit 1, "[nf-core/circrna] error: 'condition' is not a column name in the phenotype file.\n\nThe response variable must be named 'condition', please refer to the usage documentation online"

        def condition  = row.condition.matches('NA') ? 'NA' : row.condition

        if(condition == '') exit 1, "[nf-core/circrna] error: Invalid phenotype file, condition column contains empty cells."
        if(condition.matches('NA')) exit 1, "[nf-core/circrna] error: NA value in phenotype condition column."

        }
        .toList()

        return Channel.value(file(pheno))
}

/*
================================================================================
                            nf-core functions
================================================================================
*/

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/circrna v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

// handle multiqc_channels
//if(params.trim_fastq == false){
//    ch_multiqc_report = multiqc_trim_out
//}else{
//    ch_multiqc_report = multiqc_raw_out
//}

// Completion e-mail notification
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/circrna] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/circrna] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/circrna] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/circrna] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$workflow.projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$workflow.projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$workflow.projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$workflow.projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/circrna] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
            mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/circrna] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset  = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/circrna]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/circrna]${c_red} Pipeline completed with errors${c_reset}-"
    }
}



def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
