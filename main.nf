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

/*
================================================================================
                                Help Flags
================================================================================
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    ====================================================
                    circrna v${workflow.manifest.version}
   ====================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/circrna -profile <docker/singularity/institute> --input 'samples.csv' --input_type 'fastq' --module 'circrna_discovery' --tool 'circexplorer2'

   Mandatory arguments:
      -profile                        [str] Configuration profile to use. If selecting multiple, provide as a comma seperated list.
                                            Note that later profiles overwrite earlier profiles, order is important!
                                            Available: docker, singularity, <institute>

      --module                        [str] Specify analysis module(s) to run. Can use multiple (comma separated), must include 'circrna_discovery'
                                            Available: circrna_discovery, mirna_prediction, differential_expression.
                                            Default: ${params.module}

      --tool                          [str] Specify which circRNA quantification tool to use.
                                            If selecting multiple tools, provide as a comma seperated list.
                                            Available: circexplorer2, circrna_finder, ciriquant, dcc, find_circ, mapsplice.
                                            Default: ${params.tool}

      --input_type                    [str] Specify the input data type.
                                            Avaialable: fastq or bam.
                                            Default: ${params.input_type}

      --input                        [file] Either paths to FASTQ/BAM data (must be surrounded with quotes).
                                            Indicate multiple files with a suitable wildcard glob pattern i.e '*_{1,2}' for FASTQ paired end reads.

                                            OR

                                            A path to a CSV file (ending .csv) containing file URL/paths and sample IDs.
                                            Please see online documentation for a template.

      --phenotype                    [file] When differential_expression is provided to --module, a phenotype CSV file (ending .csv) must be provided.
                                            Sample IDs and response + explanatory variables from metadata must be included. Do not include irrelavant metadata.
                                            This file is used to construct the DESeq2 model design.
                                            The response variable must be named 'condition' &  wild-type/control/normal samples named 'control'.
                                            Please see online documentation for examples of a valid phenotype.csv file.

    circRNA filtering
      --bsj_reads                     [int] Define the required number of reads spanning circRNA back-splice junction for circRNAs.
                                            circRNAs with counts below [int] are discarded.
                                            Disable by setting to 0.
                                            Default: ${params.bsj_reads}

      --tool_filter                   [int] Specify the minimum number of tools circRNAs must be called by.
                                            Disable by setting to 0/1 (i.e take the union, do not apply filter).
                                            Default: ${params.tool_filter}

    Reference files
      --genome_version                [str] When running the pipeline for the first time, specify the genome version to download for the analysis.
                                            Gencode reference Fasta, GTF and annotation text files will be automatically generated.
                                            Available: 'GRCh37', 'GRCh38'.

      --fasta                        [file] Path to Gencode reference genome FASTA file. Must end with '.fa', '.fasta'.

      --gencode_gtf                  [file] Path to Gencode reference GTF file. Must end with '.gtf'.

      --gene_annotation              [file] Path to customised gene annotation file. Recommended to allow [nf-core/circrna] generate this file.

      --bowtie_index                  [dir] Path to directory containing bowtie indices.

      --bowtie2_index                 [dir] Path to directory containing bowtie2 indices.

      --bwa_index                     [dir] Path to directory containing BWA indices.

      --hisat2_index                  [dir] Path to directory containing HISAT2 indices.

      --star_index                    [dir] Path to directory containing STAR indices.

      --fasta_fai                    [file] Path to SAMtools genome index file.

    Adapter trimming
      --skip_trim                     [str] Specify whether to skip BBDUK adapter/quality trimming. Available: 'no', 'yes'

      --adapters                     [file] Path to adapters file containing sequences to trim.
                                            Requires '--k', '--ktrim' and optionally '--hdist'.
                                            Default: ${params.adapters}

      --k                             [int] Specify k-mer size to use for sequence - adapter matching.
                                            Requires '--adapters', '--ktrim' and optionally '--hdist'.
                                            Default: ${params.k}

      --ktrim                         [str] Specify which ends of reads to perform adapter trimming on.
                                            Requires '--adapters', '--k' and optionally '--hdist'.
                                            Default: ${params.ktrim}

      --hdist                         [int] Specify maximin hamming distance for reference k-mers.
                                            Default: ${params.hdist}

      --trimq                         [int] Regions within reads that fall below this average Phred score are trimmed.
                                            Requires '--qtrim'.
                                            Default: ${params.trimq}

      --qtrim                         [str] Specify which ends of reads to perform trimming on.
                                            Requires '--trimq'.
                                            Default: ${params.qtrim}

      --minlen                        [int] Filter to remove trimmed reads with length below this value.
                                            Default: ${params.minlen}

    STAR alignment
      --alignIntronMax                [int] Maximum length STAR considers for introns.
                                            Default: ${params.alignIntronMax}

      --alignIntronMin                [int] Minimum length STAR considers for introns (otherwise considered a deletion).
                                            Default: ${params.alignIntronMin}

      --alignMatesGapMax              [int] Maximum gap STAR considers between mates.
                                            Default: ${params.alignMatesGapMax}

      --alignSJDBoverhangMin          [int] Minimum overhang for annotated junctions.
                                            Default: ${params.alignSJDBoverhangMin}

      --alignSJoverhangMin            [int] Minimum overhang for unannotated junctions.
                                            Default: ${params.alignSJoverhangMin}

      --alignSoftClipAtReferenceEnds  [str] Allow soft-clipping of alignments past the ends of chromosomes.
                                            Default: ${params.alignSoftClipAtReferenceEnds}

      --alignTranscriptsPerReadNmax   [int] Max number of alignments per read to consider.
                                            Default: ${params.alignTranscriptsPerReadNmax}

      --chimJunctionOverhangMin       [int] Minimum overhang for a chimeric junction.
                                            Default: ${params.chimJunctionOverhangMin}

      --chimScoreMin                  [int] Minimum total score (summed) of chimeric segments.
                                            Default: ${params.chimScoreMin}

      --chimScoreSeparation           [int] Minimum difference between the best chimeric score and the next one.
                                            Default: ${params.chimScoreSeparation}

      --genomeLoad                    [str] Mode of shared memory usage for genome files.
                                            Available: 'LoadAndKeep', 'LoadAndRemove', 'LoadAndExit', 'Remove', 'NoSharedMemory'.
                                            Default: ${params.genomeLoad}

      --limitSjdbInsertNsj            [int] Maximum number of junctions to be inserted on the fly during mapping stage.
                                            Default: ${params.limitSjdbInsertNsj}

      --outFilterMatchNminOverLread [float] Output alignment if ratio of matched bases relative to read length is >= value.
                                            Default: ${params.outFilterMatchNminOverLread}

      --outFilterMismatchNoverLmax  [float] Output alignment if ratio of mismatched based relative to mapped read length is <= value.
                                            Default: ${params.outFilterMismatchNoverLmax}

      --outFilterMultimapNmax         [int] Maximum number of multiple alignments permitted for read.
                                            Default: ${params.outFilterMultimapNmax}

      --outFilterMultimapScoreRange   [int] Score range below the maximum score for multimapping alignments.
                                            Default: ${params.outFilterMultimapScoreRange}

      --outSJfilterOverhangMin        [str] Minimum overhang length for novel splice junctions. 4 integers provided in string.
                                            From left to right, integers sepcify minimum overhang length for splice junction sites:
                                            1. non-canonical motifs
                                            2. GT/AG and CT/AC motifs
                                            3. GC/AG and CT/GC motifs
                                            4. AT/AC and GT/AT motifs
                                            Default: ${params.outSJfilterOverhangMin}

      --sjdbOverhang                  [int] Specify length of donor/acceptor sequence flanking junctions (Index generation step).
                                            Default: ${params.sjdbOverhang}

      --sjdbScore                     [int] Alignment score for alignments that traverse database junctions.
                                            Default: ${params.sjdbScore}

      --winAnchorMultimapNmax         [int] Maximum number of loci anchors are allowed map to.
                                            Default: ${params.winAnchorMultimapNmax}

    Other options:
      --outdir                        [dir] The output directory where the results will be saved.
                                            Default: ${params.outdir}
      --publish_dir_mode             [list] Mode for publishing results in the output directory (only one)
                                            Available: symlink, rellink, link, copy, copyNoFollow, move
                                            Default: copy

   For a full description of the parameters, visit [nf-core/circrna] homepage (https://nf-co.re/circrna).
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Small console separator to make it easier to read errors after launch
println ""

/*
================================================================================
                          Check parameters
================================================================================
*/

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/circrna --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

// Check Tools selected
toolList = defineToolList()
tool = params.tool ? params.tool.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tool, toolList)) exit 1, "[nf-core/circrna] error: Unknown tool selected, see --help for more information."

// Check Modules selected
moduleList = defineModuleList()
module = params.module ? params.module.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(module, moduleList)) exit 1, "[nf-core/circrna] error: Unknown module selected, see --help for more information."

// Check outdir not empty string.
if(params.outdir == ''){
   exit 1, "[nf-core/circrna] error: --outdir was not supplied, please provide a output directory to publish workflow results."
}

// Check input type not empty.
if(params.input_type == ''){
   exit 1, "[nf-core/circrna] error: --input_type was not supplied, please select 'fastq' or 'bam'."
}

// Check Genome version
if(params.genome_version == ''){
   exit 1, "[nf-core/circrna] error: --genome_version was not supplied, please select 'GRCh37' or 'GRCh38'."
}

// Check proper versions supplied
if(params.genome_version){

   GenomeVersions = defineGenomeVersions()

   Channel
         .value(params.genome_version)
         .map{ it ->

               if(!GenomeVersions.contains(it)){
                  exit 1, "[nf-core/circrna] error: Incorrect genome version (${params.genome_version}) supplied.\n\nPlease select 'GRCh37' or 'GRCh38'"
               }
         }
}

/*
 * The below Reference Genome / index parameters are allowed to be empty (they will be generated if empty)
 * Mainly concerned about valid file extensions when provided.
 */

// Check Fasta
if(params.fasta && !has_extension(params.fasta, ".fa") && !has_extension(params.fasta, ".fasta")){
   exit 1, "[nf-core/circrna] error: Reference Fasta file provided (${params.fasta}) is not valid, Fasta file should have the extension '.fa' or '.fasta'."
}

// Check GTF
if(params.gencode_gtf && !has_extension(params.gencode_gtf, ".gtf")){
   exit 1, "[nf-core/circrna] error: Reference GTF file provided (${params.gencode_gtf}) is not valid, GTF file should have the extension '.gtf'."
}

// Check Fasta fai
if(params.fasta_fai && !has_extension(params.fasta_fai, ".fai")){
   exit 1, "[nf-core/circrna] error: Fasta index file provided (${params.fasta_fai}) is not valid, Fasta index files should have the extension '.fai'."
}

// Check BWA index
if(params.bwa_index){

   bwa_path_files = params.bwa_index + "/*"

   Channel
         .fromPath(bwa_path_files, checkIfExists: true)
         .flatten()
         .map{ it ->

               if(!has_extension(it, ".ann") && !has_extension(it, ".amb") && !has_extension(it, ".bwt") && !has_extension(it, ".pac") && !has_extension(it, ".sa")){
                  exit 1, "[nf-core/circrna] error: BWA index file ($it) has an incorrect extension. Are you sure they are BWA indices?"
               }
         }
}

// Check Bowtie index
if(params.bowtie_index){

   bowtie_path_files = params.bowtie_index + "/*"

   Channel
         .fromPath(bowtie_path_files, checkIfExists: true)
         .flatten().view()
         .map{ it ->

               if(!has_extension(it, ".ebwt")){
                  exit 1, "[nf-core/circrna] error: Bowtie index file ($it) has an incorrect extension. Are you sure they are Bowtie(1) indices?"
               }
         }
}

// Check Bowtie 2 index
if(params.bowtie2_index){

   bowtie2_path_files = params.bowtie2_index + "/*"

   Channel
         .fromPath(bowtie2_path_files, checkIfExists: true)
         .flatten()
         .map{ it ->

               if(!has_extension(it, ".bt2")){
                  exit 1, "[nf-core/circrna] error: Bowtie 2 index file ($it) has an incorrect extension. Are you sure they are Bowtie 2 indices?"
               }
          }
}


// Check HISAT2 index
if(params.hisat2_index){

   hisat2_path_files = params.hisat2_index + "/*"

   Channel
         .fromPath(hisat2_path_files, checkIfExists: true)
         .flatten()
         .map{ it ->

               if(!has_extension(it, ".ht2")){
                  exit 1, "[nf-core/circrna] error: HISAT2 index file ($it) has an incorrect extension. Are you sure they are HISAT2 indices?"
               }
          }
}

// Check STAR index
if(params.star_index){

   starList = defineStarFiles()

   star_path_files = params.star_index + "/*"

   Channel
         .fromPath(star_path_files, checkIfExists: true)
         .flatten()
         .map{ it -> it.getName()}
         .collect()
         .flatten()
         .map{ it ->

               if(!starList.contains(it)){
                  exit 1, "[nf-core/circrna] error: Incorrect index file ($it) is in the STAR directory provided.\n\nPlease check your STAR indices are valid:\n$starList."
               }
          }
}

// Check phenotype file

// Check it is a valid file & stage path:
pheno_path = null
if(params.phenotype && (has_extension(params.phenotype, ".csv"))){
   pheno_path = params.phenotype
}else{
   exit 1, "[nf-core/circrna] error: Input phenotype file (${params.phenotype}) is incorrect.\n\nMust be a '.csv' file and be comma delimited. See online documentation for description + examples."
}

// Check 'condition' is a column name.
ch_phenotype = Channel.empty()
if(pheno_path){

   pheno_file = file(pheno_path)

   if(pheno_file instanceof List) exit 1, "[nf-core/circrna] error: multiple files passed to --phenotype parameter."
   if(!pheno_file.exists()) exit 1, "[nf-core/circrna] error: input phenotype file could not be found using the path provided: ${params.phenotype}"

   ch_phenotype = examine_phenotype(pheno_file)

}

// Check BBDUK params

// Check adapters
if(params.skip_trim == 'no'){
   if(!has_extension(params.adapters, ".fa") && !has_extension(params.adapters, ".fasta")){
      exit 1, "[nf-core/circrna] error: Adapters file provied (${params.adapters}) is not valid, Fasta files must have the extension '*.fa' or '*.fasta'."
   }
   if(params.adapters){
      adapters = file(params.adapters, checkIfExists: true)
   }
}

// Check all adapter trimming flags are provided
if(params.skip_trim == 'no' && params.adapters && (!params.k && !params.ktrim || !params.k && params.ktrim || params.k && !params.ktrim)){
   exit 1, "[nf-core/circrna] error: Adapter file provided for trimming but missing values for '--k' and/or '--ktrim'.Please provide values for '--k' and '--ktrim'.\n\nPlease check the parameter documentation online."
}

// Check all quality trimming flags are provided
if(params.skip_trim == 'no' && (params.trimq && !params.qtrim || !params.trimq && params.qtrim)){
   exit 1, "[nf-core/circrna] error: Both '--trimq' and '--qtrim' are required to perform quality filtering - only one has been provided.\n\nPlease check the parameter documentation online."
}

// Check filtering params
tools_selected = tool.size()

// Check it is an integer
if(tools_selected > 1 && !isValidInteger(params.tool_filter)){
  exit 1, "[nf-core/circrna] error: The parameter '--tool_filter' ($params.tool_filter) is not an integer. Do not wrap in quotes.\n\nPlease check the documentation online."
}

// Check it does not exceed number of tools selected
if(tools_selected > 1 && params.tool_filter > tools_selected){
  exit 1, "[nf-core/circrna] error: The parameter '--tool_filter' (${params.tool_filter}) exceeds the number of tools selected (${params.tool}). Please select a value less than or equal to the number of quantification tools selected ($tools_selected).\n\nPlease check the help documentation."
}

// Check BSJ reads param
if(!isValidInteger(params.bsj_reads)){
exit 1, "[nf-core/circrna] error: The parameter '--bsj_reads' ($params.bsj_reads) is not an integer. Do not wrap in quotes.\n\nPlease check the documentation online."
}


/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName

log.info nfcoreHeader()
def summary = [:]
if (workflow.revision)        summary['Pipeline Release']    = workflow.revision
summary['Run Name']          = custom_runName ?: workflow.runName
if (workflow.containerEngine) summary['Container'] = "${workflow.containerEngine} - ${workflow.container}"
summary['Max Resources']     = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
summary['Config Files']   = workflow.configFiles.join(', ')
summary['Launch dir']  = workflow.launchDir
summary['Output dir']  = params.outdir
summary['Publish dir mode']  = params.publish_dir_mode
summary['Working dir'] = workflow.workDir
summary['Script dir']  = workflow.projectDir
summary['User']        = workflow.userName

summary['Input']             = params.input
summary['Input type']        = params.input_type
summary['circRNA tool(s)']   = params.tool
summary['modules']           = params.module
if('differential_expression' in module) summary['Phenotype design'] = params.phenotype
summary['BSJ filter']        = params.bsj_reads
if(tools_selected > 1) summary['Tool filter'] = params.tool_filter

summary['Genome version'] = params.genome_version
if(params.fasta)           summary['Reference FASTA']   = params.fasta
if(params.gencode_gtf)     summary['Reference GTF']     = params.gencode_gtf
if(params.gene_annotation) summary['Custom annotation'] = params.gene_annotation
if(params.bowtie_index)    summary['Bowtie indices']    = params.bowtie_index
if(params.bowtie2_index)   summary['Bowtie2 indices']   = params.bowtie2_index
if(params.bwa_index)       summary['BWA indices']       = params.bwa_index
if(params.fasta_fai)       summary['SAMtools index']    = params.fasta_fai
if(params.hisat2_index)    summary['HISAT2 indices']    = params.hisat2_index
if(params.star_index)      summary ['STAR indices']     = params.star_index

summary['Skip BBDUK']     = params.skip_trim
if(params.skip_trim == 'no'){
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

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

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
                          Export software versions
================================================================================
*/

process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: {it.indexOf(".csv") > 0 ? it : null}

    output:
        file 'software_versions_mqc.yaml' into ch_software_versions_yaml
        file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    bbduk.sh | grep 'Last modified' | cut -d' ' -f 3-99 &> v_bbduk.txt 2>&1 || true
    bedtools --version &> v_bedtools.txt 2>&1 || true
    bowtie --version | awk -v OFS=' ' '{print \$3}' | head -n 1 &> v_bowtie.txt 2>&1 || true
    bowtie2 --version | awk -v OFS=' ' '{print \$3}' | head -n 1 &> v_bowtie2.txt 2>&1 || true
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' &> v_bwa.txt 2>&1 || true
    CIRCexplorer2 --version &> v_circexplorer2.txt 2>&1 || true
    CIRIquant --version &> v_ciriquant.txt 2>&1 || true
    hisat2 --version &> v_hisat2.txt 2>&1 || true
    java --version &> v_java.txt 2>&1 || true
    mapsplice.py --version &> v_mapsplice.txt 2>&1 || true
    miranda -v | grep 'Algorithm' | cut -d' ' -f 1,2 &> v_miranda.txt 2>&1 || true
    perl --version | grep 'x86' | cut -d' ' -f1-9 &> v_perl.txt 2>&1 || true
    picard SamToFastq --version &> v_picard.txt 2>&1 || true
    pip --version | cut -d' ' -f 1,2,5,6 &> v_pip.txt 2>&1 || true
    python --version &> v_python.txt 2>&1 || true
    R --version | head -n 1 &> v_R.txt || true
    RNAfold --version | cut -d' ' -f2 &> v_viennarna.txt || true
    samtools --version &> v_samtools.txt 2>&1 || true
    STAR --version &> v_star.txt 2>&1 || true
    stringtie --version &> v_stringtie.txt 2>&1 || true
    targetscan_70.pl --version &> v_targetscan.txt 2>&1 || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
================================================================================
                          Begin Workflow
================================================================================
*/

/*
================================================================================
                          Download Files
================================================================================
*/

process download_fasta{

    publishDir "${params.outdir}/circrna_discovery/reference", mode: params.publish_dir_mode

    output:
        file("*.fa") into fasta_downloaded

    when: !params.fasta

    shell:
    if(params.genome_version == 'GRCh37'){
       $/
       wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
       gunzip GRCh37.primary_assembly.genome.fa.gz
       mv GRCh37.primary_assembly.genome.fa GRCh37.fa.tmp
       sed 's/\s.*$//' GRCh37.fa.tmp > GRCh37.fa
       rm GRCh37.fa.tmp
       /$
    }else if(params.genome_version == 'GRCh38'){
       $/
       wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
       gunzip GRCh38.primary_assembly.genome.fa.gz
       mv GRCh38.primary_assembly.genome.fa GRCh38.fa.tmp
       sed 's/\s.*$//' GRCh38.fa.tmp > GRCh38.fa
       rm GRCh38.fa.tmp
       /$
    }
}

ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded

process download_gtf{

    publishDir "${params.outdir}/circrna_discovery/reference", mode: params.publish_dir_mode

    output:
        file("*.gtf") into gtf_downloaded

    when: !params.gencode_gtf

    shell:
    if(params.genome_version == 'GRCh37'){
       $/
       wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz
       gunzip gencode.v34lift37.annotation.gtf.gz
       mv gencode.v34lift37.annotation.gtf GRCh37.gtf
       /$
    }else if(params.genome_version == 'GRCh38'){
       $/
       wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
       gunzip gencode.v34.primary_assembly.annotation.gtf.gz
       mv gencode.v34.primary_assembly.annotation.gtf GRCh38.gtf
       /$
    }
}

ch_gencode_gtf = params.gencode_gtf ? Channel.value(file(params.gencode_gtf)) : gtf_downloaded

process create_gene_annotation{

    publishDir "${params.outdir}/circrna_discovery/reference/", mode: params.publish_dir_mode

    input:
        file(gtf) from ch_gencode_gtf

    output:
        file("*.txt") into gene_annotation_created

    when: !params.gene_annotation

    shell:
    $/
    gtfToGenePred -genePredExt -geneNameAsName2 !{gtf} !{params.genome_version}.genepred
    perl -alne '$"="\t";print "@F[11,0..9]"' !{params.genome_version}.genepred > !{params.genome_version}.txt
    /$
}

ch_gene_annotation = params.gene_annotation ? Channel.value(file(params.gene_annotation)) : gene_annotation_created

process download_mirbase{

    errorStrategy 'retry'
    maxRetries 10

    publishDir "${params.outdir}/mirna_prediction/assets",mode: params.publish_dir_mode

    output:
        file("hsa_mature.fa") into miranda_miRs

    when: 'mirna_prediction' in module

    script:
    """
    wget --no-check-certificate ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
    gunzip mature.fa.gz
    grep "sapiens" -A1 mature.fa | awk '!/--/' > hsa_mature.fa
    """
}

process download_targetscan{

    errorStrategy 'retry'
    maxRetries 10

    publishDir "${params.outdir}/mirna_prediction/assets", mode: params.publish_dir_mode

    output:
        file("hsa_miR.txt") into targetscan_miRs
        file("hsa_miR_for_context_scores.txt") into targetscan_miRs_context_scores

    when: 'mirna_prediction' in module

    script:
    """
    wget --no-check-certificate http://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip
    jar xvf miR_Family_Info.txt.zip
    grep 9606 miR_Family_Info.txt > hsa_miR_Family_Info.txt
    awk -v OFS="\t" '{print \$1, \$2, \$3}' hsa_miR_Family_Info.txt > hsa_miR.txt
    awk -v OFS="\t" '{print \$1, \$3, \$4, \$5}' hsa_miR.txt > hsa_miR_for_context_scores.txt
    """
}


/*
================================================================================
                          Create Genome Indices
================================================================================
*/

process samtools_index{

    publishDir "${params.outdir}/circrna_discovery/index/samtools", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.fai") into fasta_fai_built

    when: !(params.fasta_fai)

    script:
    """
    samtools faidx $fasta
    """
}

ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fasta_fai_built

process bwa_index{

    publishDir "${params.outdir}/circrna_discovery/index/bwa", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta.baseName}.*") into bwa_built
        val("$launchDir/${params.outdir}/circrna_discovery/index/bwa") into bwa_path

    when: !(params.bwa_index) && 'ciriquant' in tool && 'circrna_discovery' in module

    script:
    """
    bwa index -a bwtsw $fasta -p ${fasta.baseName}
    """
}

ch_bwa_index = params.bwa_index ? Channel.value(params.bwa_index) : bwa_path

process hisat2_index{

    publishDir "${params.outdir}/circrna_discovery/index/hisat2", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta.baseName}.*.ht2") into hisat2_built
        val("$launchDir/${params.outdir}/circrna_discovery/index/hisat2") into hisat2_path

    when: !(params.hisat2_index) && ( 'ciriquant' in tool || 'differential_expression' in module )

    script:
    """
    hisat2-build $fasta ${fasta.baseName}
    """
}

ch_hisat2_index = params.hisat2_index ? Channel.value(params.hisat2_index) : hisat2_path

process star_index{

    publishDir "${params.outdir}/circrna_discovery/index", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta
        file(gtf) from ch_gencode_gtf

    output:
        file("STAR") into star_built

    when: !(params.star_index) && ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --runThreadN ${task.cpus} \
    --sjdbOverhang ${params.sjdbOverhang} \
    --sjdbGTFfile $gtf \
    --genomeDir STAR/ \
    --genomeFastaFiles $fasta
    """
}

ch_star_index = params.star_index ? Channel.value(file(params.star_index)) : star_built

process bowtie_index{

    publishDir "${params.outdir}/circrna_discovery/index/bowtie", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta

    output:
        file ("${fasta.baseName}.*") into bowtie_built

    when: !(params.bowtie_index) && ('mapsplice' in tool || 'uroborus' in tool) && 'circrna_discovery' in module

    script:
    """
    bowtie-build $fasta ${fasta.baseName}
    """
}

bowtie_path_files = params.bowtie_index + "/*"
ch_bowtie_index = params.bowtie_index ? Channel.value(file(bowtie_path_files)) : bowtie_built

process bowtie2_index{

    publishDir "${params.outdir}/circrna_discovery/index/bowtie2", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta

    output:
        file ("${fasta.baseName}.*") into bowtie2_built

    when: !(params.bowtie2_index) && ('find_circ' in tool || 'uroborus' in tool) && 'circrna_discovery' in module

    script:
    """
    bowtie2-build $fasta ${fasta.baseName}
    """
}

bowtie2_path_files = params.bowtie2_index + "/*"
ch_bowtie2_index = params.bowtie2_index ? Channel.value(file(bowtie2_path_files)) : bowtie2_built

/*
================================================================================
                       Misc. circRNA requirements
================================================================================
*/


process split_fasta{

    publishDir "${params.outdir}/circrna_discovery/reference/chromosomes", mode: params.publish_dir_mode

    input:
        file(fasta) from ch_fasta

    output:
        path("*.fa", includeInputs:true) into split_fasta
        val("${launchDir}/${params.outdir}/circrna_discovery/reference/chromosomes") into ch_fasta_chr

    when: ('mapsplice' in tool || 'find_circ' in tool) && 'circrna_discovery' in module

    shell:
    '''
    ## Add catch for test data (uses only 1 chr, no action needed)
    n_chr=$(grep '>' !{fasta} | wc -l)

    if [[ $n_chr -gt 1 ]];
    then
      awk '/^>/ {F=substr($0, 2, length($0))".fa"; print >F;next;} {print >> F;}' < !{fasta}
      rm !{fasta}
    else
      :
    fi
    '''
}

process ciriquant_yml{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/ciriquant", mode: params.publish_dir_mode

    input:
        file(gencode_gtf) from ch_gencode_gtf
        file(fasta) from ch_fasta
        val(bwa_path) from ch_bwa_index
        val(hisat2_path) from ch_hisat2_index

    output:
        file("travis.yml") into ch_ciriquant_yml

    when: 'ciriquant' in tool && 'circrna_discovery' in module

    script:
    index_prefix = fasta.toString() - ~/.fa/
    fasta_path = fasta.toRealPath()
    gencode_gtf_path = gencode_gtf.toRealPath()
    """
    export bwa=`whereis bwa | cut -f2 -d':'`
    export hisat2=`whereis hisat2 | cut -f2 -d':'`
    export stringtie=`whereis stringtie | cut -f2 -d':'`
    export samtools=`whereis samtools | cut -f2 -d':' | awk '{print \$1}'`

    touch travis.yml
    printf "name: ciriquant\n\
    tools:\n\
     bwa: \$bwa\n\
     hisat2: \$hisat2\n\
     stringtie: \$stringtie\n\
     samtools: \$samtools\n\n\
    reference:\n\
     fasta: ${fasta_path}\n\
     gtf: ${gencode_gtf_path}\n\
     bwa_index: ${bwa_path}/${index_prefix}\n\
     hisat_index: ${hisat2_path}/${index_prefix}" >> travis.yml
    """
}

/*
================================================================================
                          Process Input Data
================================================================================
*/

// Check inputs

if(params.input == null){
   exit 1, "[nf-core/circrna] error: --input was not supplied! Please check '--help' or documentation for details"
}

csv_path = null
if(params.input && (has_extension(params.input, "csv"))) csv_path = params.input

ch_input_sample = Channel.empty()
if(csv_path){

   csv_file = file(csv_path)
   println csv_file
   if(csv_file instanceof List) exit 1, "[nf-core/circrna] error: can only accept one CSV file per run."
   if(!csv_file.exists()) exit 1, "[nf-core/circrna] error: input CSV file could not be found using the path provided: ${params.input}"

   ch_input_sample = extract_data(csv_path)

}else if(params.input && !has_extension(params.input, "csv")){

   log.info ""
   log.info "Input data log info:"
   log.info "No input sample CSV file provided, attempting to read from path instead."
   log.info "Reading input data from path: ${params.input}\n"
   log.info ""
   inputSample = retrieve_input_paths(params.input, params.input_type)
   ch_input_sample = inputSample

}else exit 1, "[nf-core/circrna] error: --input file(s) not correctly supplied or improperly defined, see '--help' flag or documentation"

// Process bam file / stage fastq

if(params.input_type == 'bam'){

   process bam_to_fq{

        publishDir "${params.outdir}/quality_control/preprocessing/bamtofastq", mode: params.publish_dir_mode

        input:
            tuple val(base), file(bam) from ch_input_sample

        output:
            tuple val(base), file('*.fq.gz') into fastq_built

        script:
        """
        picard "-Xmx${task.memory.toGiga()}g" \
        SamToFastq \
        I=$bam \
        F=${base}_R1.fq.gz \
        F2=${base}_R2.fq.gz \
        VALIDATION_STRINGENCY=LENIENT
        """
   }

   (fastqc_reads, trimming_reads, raw_reads) = fastq_built.into(3)

}else if(params.input_type == 'fastq'){

   (fastqc_reads, trimming_reads, raw_reads) = ch_input_sample.into(3)

}else exit 1, "[nf-core/circrna] error: --input_type must be one of 'fastq' or 'bam'."

// FASTQC on raw data

process FastQC {

    label 'py3'

    publishDir "${params.outdir}/quality_control/fastqc/raw", mode: params.publish_dir_mode

    input:
        tuple val(base), file(fastq) from fastqc_reads

    output:
        file("*.{html,zip}") into fastqc_raw

    script:
    """
    fastqc -q $fastq
    """
}

// MultiQC of the Raw Data

// collect software_versions, workflow summary
(software_versions_raw, software_versions_trim) = ch_software_versions_yaml.into(2)
(workflow_summary_raw, workflow_summary_trim) = ch_workflow_summary.into(2)
process multiqc_raw {

    publishDir "${params.outdir}/quality_control/multiqc", mode: params.publish_dir_mode

    label 'py3'

    input:
        file(htmls) from fastqc_raw.collect()
        file ('software_versions/*') from software_versions_raw.collect()
        file workflow_summary from workflow_summary_raw.collectFile(name: "workflow_summary_mqc.yaml")

    output:
        file("Raw_Reads_MultiQC.html") into multiqc_raw_out

    script:
    """
    multiqc -i "Raw_Reads_MultiQC" -b "nf-circ pipeline" -n "Raw_Reads_MultiQC.html" .
    """
}

// BBDUK

if(params.skip_trim == 'no'){

   process bbduk {

       publishDir "${params.outdir}/quality_control/preprocessing/BBDUK", pattern: "*fq.gz", mode: params.publish_dir_mode

       input:
           tuple val(base), file(fastq) from trimming_reads
           path adapters from params.adapters

       output:
           tuple val(base), file('*.trim.fq.gz') into trim_reads_ch
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
       #uncoment for nf-core, NUIG slurm does not have memory param
       #bbduk.sh -Xmx${task.memory.toGiga()}g \
       bbduk.sh -Xmx4g \
       threads=${task.cpus} \
       in1=${fastq[0]} \
       in2=${fastq[1]} \
       out1=${base}_R1.trim.fq.gz \
       out2=${base}_R2.trim.fq.gz \
       $adapter \
       $k \
       $ktrim \
       $trimq \
       $qtrim \
       $minlen \
       stats=${base}_BBDUK.txt
       """
   }

   // trimmed reads into 2 channels:
   (fastqc_trim_reads, aligner_reads) = trim_reads_ch.into(2)

   process FastQC_trim {

       label 'py3'

       publishDir "${params.outdir}/quality_control/fastqc/trimmed", mode: params.publish_dir_mode

       input:
           tuple val(base), file(fastq) from fastqc_trim_reads

       output:
           file ("*.{html,zip}") into fastqc_trimmed

       script:
       """
       fastqc -q $fastq
       """
   }

   process multiqc_trim {

       publishDir "${params.outdir}/quality_control/multiqc", mode: params.publish_dir_mode

       label 'py3'

       input:
           file(htmls) from fastqc_trimmed.collect()
           file(bbduk_stats) from bbduk_stats_ch.collect()
           file ('software_versions/*') from software_versions_trim.collect()
           file workflow_summary from workflow_summary_trim.collectFile(name: "workflow_summary_mqc.yaml")

       output:
           file("Trimmed_Reads_MultiQC.html") into multiqc_trim_out

       script:
       """
       multiqc -i "Trimmed_Reads_MultiQC" -b "nf-circ pipeline" -n "Trimmed_Reads_MultiQC.html" .
       """
   }
}else if(params.skip_trim == 'yes'){
   aligner_reads = raw_reads
}

// Stage Aligner read channels
(star_pass1_reads, star_pass2_reads, find_circ_reads, ciriquant_reads, mapsplice_reads, uroborus_reads, dcc_mate1_reads, dcc_mate2_reads, hisat2_reads) = aligner_reads.into(9)

/*
================================================================================
                             circRNA Discovery
================================================================================
*/

// STAR 1st Pass

process STAR_1PASS{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/STAR/1st_Pass", pattern: "${base}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(reads) from star_pass1_reads
        val(star_idx) from ch_star_index

    output:
        file("${base}/*SJ.out.tab") into sjdb_ch
        file("${base}") into star_1st_pass_output

    when: ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    STAR \
    --alignIntronMax ${params.alignIntronMax} \
    --alignIntronMin ${params.alignIntronMin} \
    --alignMatesGapMax ${params.alignMatesGapMax} \
    --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \
    --alignSJoverhangMin ${params.alignSJoverhangMin} \
    --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \
    --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \
    --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \
    --chimOutType Junctions SeparateSAMold \
    --chimScoreMin ${params.chimScoreMin} \
    --chimScoreSeparation ${params.chimScoreSeparation} \
    --chimSegmentMin ${params.chimSegmentMin} \
    --genomeDir ${star_idx} \
    --genomeLoad ${params.genomeLoad} \
    --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \
    --outFileNamePrefix ${base}/${base}. \
    --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \
    --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \
    --outFilterMultimapNmax ${params.outFilterMultimapNmax} \
    --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \
    --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \
    --outFilterType BySJout \
    --outReadsUnmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \
    ${readFilesCommand} \
    --readFilesIn ${reads} \
    --runThreadN ${task.cpus} \
    --sjdbScore ${params.sjdbScore} \
    --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

process sjdbFile{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/STAR/SJFile", pattern: "*SJFile.tab", mode: params.publish_dir_mode

    input:
        file(sjdb) from sjdb_ch

    output:
        file("*SJFile.tab") into sjdbfile_ch

    when: ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    shell:
    '''
    base=$(basename !{sjdb} .SJ.out.tab)
    awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' !{sjdb} > ${base}.SJFile.tab
    '''
}

(sjdbfile_pass2, sjdbfile_mate1, sjdbfile_mate2) = sjdbfile_ch.into(3)

process STAR_2PASS{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/STAR/2nd_Pass", pattern: "${base}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(reads) from star_pass2_reads
        file(sjdbfile) from sjdbfile_pass2.collect()
        val(star_idx) from ch_star_index

    output:
        tuple val(base), file("${base}/${base}.Chimeric.out.junction") into circexplorer2_input
        tuple val(base), file("${base}") into circrna_finder_input, dcc_pairs

    when: ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    STAR \
    --alignIntronMax ${params.alignIntronMax} \
    --alignIntronMin ${params.alignIntronMin} \
    --alignMatesGapMax ${params.alignMatesGapMax} \
    --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \
    --alignSJoverhangMin ${params.alignSJoverhangMin} \
    --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \
    --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \
    --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \
    --chimOutType Junctions SeparateSAMold \
    --chimScoreMin ${params.chimScoreMin} \
    --chimScoreSeparation ${params.chimScoreSeparation} \
    --chimSegmentMin ${params.chimSegmentMin} \
    --genomeDir ${star_idx} \
    --genomeLoad ${params.genomeLoad} \
    --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \
    --outFileNamePrefix ${base}/${base}. \
    --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \
    --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \
    --outFilterMultimapNmax ${params.outFilterMultimapNmax} \
    --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \
    --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \
    --outFilterType BySJout \
    --outReadsUnmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \
    ${readFilesCommand} \
    --readFilesIn ${reads} \
    --runThreadN ${task.cpus} \
    --sjdbFileChrStartEnd ${sjdbfile} \
    --sjdbScore ${params.sjdbScore} \
    --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

// CIRCexplorer2

process circexplorer2_star{

    publishDir "${params.outdir}/circrna_discovery/filtered_outputs/circexplorer2", pattern: "*_circexplorer2.bed", mode: params.publish_dir_mode
    publishDir "${params.outdir}/circrna_discovery/tool_outputs/circexplorer2", pattern: "${base}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(chimeric_reads) from circexplorer2_input
        file(fasta) from ch_fasta
        file(gene_annotation) from ch_gene_annotation

    output:
        tuple val(base), file("${base}_circexplorer2.bed") into circexplorer2_results
        tuple val(base), file("${base}") into circexplorer2_raw


    when: 'circexplorer2' in tool && 'circrna_discovery' in module

    script:
    """
    mkdir -p ${base}

    CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}/${base}.STAR.junction.bed

    CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}/${base}.STAR.junction.bed -o ${base}/${base}.txt

    awk '{if(\$13 >= ${params.bsj_reads}) print \$0}' ${base}/${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_circexplorer2.bed
    """
}

// circRNA_finder

process circrna_finder{

    publishDir "${params.outdir}/circrna_discovery/filtered_outputs/circrna_finder", pattern: '*_circrna_finder.bed', mode: params.publish_dir_mode
    publishDir "${params.outdir}/circrna_discovery/tool_outputs/circrna_finder/${base}", pattern: "{*filteredJunctions*,*.Chimeric.out.sorted.*}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(star_dir) from circrna_finder_input

    output:
        tuple val(base), file("${base}_circrna_finder.bed") into circrna_finder_results
        tuple val(base), file("{*filteredJunctions*,*.Chimeric.out.sorted.*}") into circrna_finder_raw

    when: 'circrna_finder' in tool && 'circrna_discovery' in module

    script:
    """
    postProcessStarAlignment.pl --starDir ${star_dir}/ --outDir ./

    awk '{if(\$5 >= ${params.bsj_reads}) print \$0}' ${base}.filteredJunctions.bed | awk  -v OFS="\t" -F"\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_circrna_finder.bed
    """
}

// DCC

process dcc_mate1{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/dcc/${base}", pattern: "mate1", mode: params.publish_dir_mode

    input:
        tuple val(base), file(reads) from dcc_mate1_reads
        file(sjdbfile) from sjdbfile_mate1.collect()
        val(star_idx) from ch_star_index

    output:
        tuple val(base), file("mate1") into dcc_mate1

    when: 'dcc' in tool && 'circrna_discovery' in module

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    STAR \
    --alignIntronMax ${params.alignIntronMax} \
    --alignIntronMin ${params.alignIntronMin} \
    --alignMatesGapMax ${params.alignMatesGapMax} \
    --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \
    --alignSJoverhangMin ${params.alignSJoverhangMin} \
    --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \
    --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \
    --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \
    --chimOutType Junctions SeparateSAMold \
    --chimScoreMin ${params.chimScoreMin} \
    --chimScoreSeparation ${params.chimScoreSeparation} \
    --chimSegmentMin ${params.chimSegmentMin} \
    --genomeDir ${star_idx} \
    --genomeLoad ${params.genomeLoad} \
    --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \
    --outFileNamePrefix mate1/${base}. \
    --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \
    --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \
    --outFilterMultimapNmax ${params.outFilterMultimapNmax} \
    --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \
    --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \
    --outFilterType BySJout \
    --outReadsUnmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \
    ${readFilesCommand} \
    --readFilesIn ${reads} \
    --runThreadN ${task.cpus} \
    --sjdbFileChrStartEnd ${sjdbfile} \
    --sjdbScore ${params.sjdbScore} \
    --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

process dcc_mate2{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/dcc/${base}", pattern: "mate2", mode: params.publish_dir_mode

    input:
        tuple val(base), file(reads) from dcc_mate2_reads
        file(sjdbfile) from sjdbfile_mate2.collect()
        val(star_idx) from ch_star_index

    output:
        tuple val(base), file("mate2") into dcc_mate2

    when: 'dcc' in tool && 'circrna_discovery' in module

	  script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    STAR \
    --alignIntronMax ${params.alignIntronMax} \
    --alignIntronMin ${params.alignIntronMin} \
    --alignMatesGapMax ${params.alignMatesGapMax} \
    --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \
    --alignSJoverhangMin ${params.alignSJoverhangMin} \
    --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \
    --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \
    --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \
    --chimOutType Junctions SeparateSAMold \
    --chimScoreMin ${params.chimScoreMin} \
    --chimScoreSeparation ${params.chimScoreSeparation} \
    --chimSegmentMin ${params.chimSegmentMin} \
    --genomeDir ${star_idx} \
    --genomeLoad ${params.genomeLoad} \
    --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \
    --outFileNamePrefix mate2/${base}. \
    --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \
    --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \
    --outFilterMultimapNmax ${params.outFilterMultimapNmax} \
    --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \
    --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \
    --outFilterType BySJout \
    --outReadsUnmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \
    ${readFilesCommand} \
    --readFilesIn ${reads} \
    --runThreadN ${task.cpus} \
    --sjdbFileChrStartEnd ${sjdbfile} \
    --sjdbScore ${params.sjdbScore} \
    --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

ch_dcc_dirs = dcc_pairs.join(dcc_mate1).join(dcc_mate2)

process dcc{

    label 'py3'

    publishDir "${params.outdir}/circrna_discovery/filtered_outputs/dcc", pattern: "${base}_dcc.bed", mode: params.publish_dir_mode
    publishDir "${params.outdir}/circrna_discovery/tool_outputs/dcc/${base}", pattern: "{*.log,*Circ*}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(pairs), file(mate1), file(mate2) from ch_dcc_dirs
        file(gtf) from ch_gencode_gtf
        file(fasta) from ch_fasta

    output:
        tuple val(base), file("${base}_dcc.bed") into dcc_results
        tuple val(base), file("${base}{.log,*.Circ*}") into dcc_raw_results

    when: 'dcc' in tool && 'circrna_discovery' in module

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

    mv CircCoordinates ${base}.CircCoordinates
    mv CircRNACount ${base}.CircRNACount
    mv *.log ${base}.log
    """
}

// find_circ

process find_anchors{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/find_circ/${base}", pattern: "{*anchors.qfa.gz,*.bam}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(fastq) from find_circ_reads
        file(fasta) from ch_fasta
        file(bowtie2_index) from ch_bowtie2_index.collect()

    output:
        tuple val(base), file("${base}_anchors.qfa.gz") into ch_anchors
        tuple val(base), file("${base}{_anchors.qfa.gz,_unmapped.bam}") into find_circ_dir

    when: 'find_circ' in tool && 'circrna_discovery' in module

    script:
    """
    bowtie2 -p ${task.cpus} --very-sensitive --mm -D 20 --score-min=C,-15,0 \
    -x ${fasta.baseName} -q -1 ${fastq[0]} -2 ${fastq[1]} \
    | samtools view -hbuS - | samtools sort --threads ${task.cpus} -m 2G - > ${base}.bam

    samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

    unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz
    """
}

process find_circ{

    publishDir "${params.outdir}/circrna_discovery/filtered_outputs/find_circ/", pattern: '*_find_circ.bed', mode: params.publish_dir_mode
    publishDir "${params.outdir}/circrna_discovery/tool_outputs/find_circ/${base}", pattern: "*.sites.*", mode: params.publish_dir_mode

    input:
        tuple val(base), file(anchors) from ch_anchors
        file(bowtie2_index) from ch_bowtie2_index.collect()
        file(fasta) from ch_fasta
        val(fasta_chr_path) from ch_fasta_chr

    output:
        tuple val(base), file("${base}_find_circ.bed") into find_circ_results
        tuple val(base), file("*.sites.*") into find_circ_raw_results

    when: 'find_circ' in tool && 'circrna_discovery' in module

    script:
    """
    bowtie2 -p ${task.cpus} --reorder --mm -D 20 --score-min=C,-15,0 -q -x ${fasta.baseName} \
    -U $anchors | python ${projectDir}/bin/find_circ.py -G $fasta_chr_path -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

    ## filtering
    grep circ ${base}.sites.bed | grep -v chrM | python ${projectDir}/bin/sum.py -2,3 | python ${projectDir}/bin/scorethresh.py -16 1 | python ${projectDir}/bin/scorethresh.py -15 2 | python ${projectDir}/bin/scorethresh.py -14 2 | python ${projectDir}/bin/scorethresh.py 7 ${params.bsj_reads} | python ${projectDir}/bin/scorethresh.py 8,9 35 | python ${projectDir}/bin/scorethresh.py -17 100000 >> ${base}.txt

    tail -n +2 ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_find_circ.bed
    """
}

// CIRIquant

process ciriquant{

    publishDir "${params.outdir}/circrna_discovery/filtered_outputs/ciriquant", pattern: "${base}_ciriquant.bed", mode: params.publish_dir_mode
    publishDir "${params.outdir}/circrna_discovery/tool_outputs/ciriquant", pattern: "${base}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(fastq) from ciriquant_reads
        file(ciriquant_yml) from ch_ciriquant_yml

    output:
        tuple val(base), file("${base}_ciriquant.bed") into ciriquant_results
        tuple val(base), file("${base}") into ciriquant_raw_dir

    when: 'ciriquant' in tool && 'circrna_discovery' in module

    script:
    """
    CIRIquant -t ${task.cpus} \
    -1 ${fastq[0]} \
    -2 ${fastq[1]} \
    --config $ciriquant_yml \
    --no-gene \
    -o ${base} \
    -p ${base}

    cp ${base}/${base}.gtf .

    ## extract counts (convert float to int)
    grep -v "#" ${base}.gtf | awk '{print \$14}' | cut -d '.' -f1 > counts
    grep -v "#" ${base}.gtf | awk -v OFS="\t" '{print \$1,\$4,\$5,\$7}' > ${base}.tmp
    paste ${base}.tmp counts > ${base}_unfilt.bed

    ## filter
    awk '{if(\$5 >= ${params.bsj_reads}) print \$0}' ${base}_unfilt.bed > ${base}_filt.bed
    grep -v '^\$' ${base}_filt.bed > ${base}_ciriquant

    ## correct offset bp position
    awk -v OFS="\t" '{\$2-=1;print}' ${base}_ciriquant > ${base}_ciriquant.bed

    rm ${base}.gtf
    """
}


// mapsplice

process mapsplice_align{

    publishDir "${params.outdir}/circrna_discovery/tool_outputs/mapsplice", pattern: "${base}", mode: params.publish_dir_mode

    input:
        tuple val(base), file(fastq) from mapsplice_reads
        val(mapsplice_ref) from ch_fasta_chr
        file(bowtie_index) from ch_bowtie_index.collect()
        file(gtf) from ch_gencode_gtf

    output:
        tuple val(base), file("${base}/fusions_raw.txt") into mapsplice_fusion
        tuple val(base), file("${base}") into mapsplice_raw

    when: 'mapsplice' in tool && 'circrna_discovery' in module

    script:
    if(fastq[0].toString().endsWith('.gz')){
       prefix = gtf.toString() - ~/.gtf/
       strip1 = fastq[0].toString() - ~/.gz/
       strip2 = fastq[1].toString() - ~/.gz/
       """
       gzip -d --force ${fastq[0]}
       gzip -d --force ${fastq[1]}

       mapsplice.py \
       -c $mapsplice_ref \
       -x $prefix \
       -1 ${strip1} \
       -2 ${strip2} \
       -p ${task.cpus} \
       --bam \
       --seglen 25 \
       --min-intron ${params.alignIntronMin} \
       --max-intron ${params.alignIntronMax} \
       --min-map-len 40 \
       --fusion-non-canonical \
       --min-fusion-distance 200 \
       --gene-gtf $gtf \
       -o $base
       """
    }else{
       prefix = gtf.toString() - ~/.gtf/
       """
       mapsplice.py \
       -c $mapsplice_ref \
       -x $prefix \
       -1 ${fastq[0]} \
       -2 ${fastq[1]} \
       -p ${task.cpus} \
       --bam \
       --seglen 25 \
       --min-intron ${params.alignIntronMin} \
       --max-intron ${params.alignIntronMax} \
       --min-map-len 40 \
       --fusion-non-canonical \
       --min-fusion-distance 200 \
       --gene-gtf $gtf \
       -o $base
       """
    }
}


process mapsplice_parse{

    publishDir "${params.outdir}/circrna_discovery/filtered_outputs/mapsplice", pattern: "*_mapsplice.bed", mode: params.publish_dir_mode

    input:
        tuple val(base), file(raw_fusion) from mapsplice_fusion
        file(fasta) from ch_fasta
        file(gene_annotation) from ch_gene_annotation

    output:
        tuple val(base), file("${base}_mapsplice.bed") into mapsplice_results

    when: 'mapsplice' in tool && 'circrna_discovery' in module

    script:
    """
    mkdir -p ${base}

    CIRCexplorer2 parse -t MapSplice $raw_fusion -b ${base}.mapsplice.junction.bed

    CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.mapsplice.junction.bed -o ${base}.txt

    awk '{if(\$13 >= ${params.bsj_reads}) print \$0}' ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_mapsplice.bed
    """
}

// UROBORUS retry with py3 container.

process tophat_align{

    input:
        tuple val(base), file(fastq) from uroborus_reads
        file(bowtie2_index) from ch_bowtie2_index.collect()
        file(fasta) from ch_fasta

    output:
        tuple val(base), file("unmapped.bam") into tophat_unmapped_bam
        tuple val(base), file("accepted_hits.bam") into tophat_accepted_hits

    when: 'uroborus' in tool && 'circrna_discovery' in module

    script:
    """
    tophat -p ${task.cpus} -o ${base} ${fasta.baseName} ${fastq[0]} ${fastq[1]}
    mv ${base}/unmapped.bam ./
    mv ${base}/accepted_hits.bam ./
    """
}


process uroborus{

    publishDir "${params.outdir}/circrna_discovery/uroborus", mode: params.publish_dir_mode

    input:
        tuple val(base), file(unmapped_bam) from tophat_unmapped_bam
        tuple val(base), file(accepted_hits) from tophat_accepted_hits
        file(bowtie_index) from ch_bowtie_index.collect()
        file(gtf) from ch_gencode_gtf
        val(uroborus_ref) from ch_fasta_chr
        file(fasta) from ch_fasta

    output:
        file("${base}.txt") into uroborus_results

    when: 'uroborus' in tool && 'circrna_discovery' in module

    script:
    """
    samtools view $unmapped_bam > unmapped.sam

    perl /opt/conda/envs/circrna/bin/UROBORUS.pl \
    -index ${fasta.baseName} \
    -gtf $gtf \
    -fasta $uroborus_ref \
    unmapped.sam $accepted_hits &> uroborus_logs.txt

    mv circRNA_list.txt ${base}.txt
    """
}

/*
================================================================================
                            Annotate circRNAs
================================================================================
*/

if(tools_selected > 1){

   // Attempted BUG fix: remainder: true, allow empty channels (null in tuple, input.1 etc in workdir)
   // Causes WARN: input caridnality does not match (due to null items), but script works.
   // No WARN if ciriquant selected in tool.
   combined_tool = ciriquant_results.join(circexplorer2_results, remainder: true).join(dcc_results, remainder: true).join(circrna_finder_results, remainder: true).join(find_circ_results, remainder: true).join(mapsplice_results, remainder: true)

   process consolidate_algorithms{

       input:
           tuple val(base), file(ciriquant), file(circexplorer2), file(dcc), file(circrna_finder), file(find_circ), file(mapsplice) from combined_tool

       output:
           file("${base}.bed") into sample_counts

       when: 'circrna_discovery' in module

       script:
       """
       ## make tool output csv file
       files=\$(ls *.bed)

       for i in \$files; do
           printf "\$i\n" >> samples.csv
       done

       ## Add catch for empty file in tool output
       bash ${projectDir}/bin/check_empty.sh

       ## Use intersection of "n" (params.tool_filter) circRNAs called by tools
       ## remove duplicate IDs, keep highest count.
       Rscript ${projectDir}/bin/consolidate_algorithms_intersection.R samples.csv $params.tool_filter

       mv combined_counts.bed ${base}.bed
       """
   }

   process get_counts_combined{

       publishDir "${params.outdir}/circrna_discovery/count_matrix", pattern: "matrix.txt", mode: params.publish_dir_mode

       input:
           file(bed) from sample_counts.collect()

       output:
           file("circRNA_matrix.txt") into circRNA_counts
           file("matrix.txt") into matrix

       when: 'circrna_discovery' in module

       script:
       """
       python ${projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
       Rscript ${projectDir}/bin/reformat_count_matrix.R
       """
   }
}else{

   single_tool = ciriquant_results.mix(circexplorer2_results, dcc_results, circrna_finder_results, find_circ_results, mapsplice_results)

   process get_counts_single{

       publishDir "${params.outdir}/circrna_discovery/count_matrix", pattern: "matrix.txt", mode: params.publish_dir_mode

       input:
           file(bed) from single_tool.collect()
           val(tool) from params.tool

       output:
           file("circRNA_matrix.txt") into circRNA_counts
           file("matrix.txt") into matrix

       when: 'circrna_discovery' in module

       script:
       """
       for b in *.bed; do
           foo=\${b%".bed"};
           bar=\${foo%"_${tool}"};
           mv \$b \${bar}.bed
       done

       python ${projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
       Rscript ${projectDir}/bin/reformat_count_matrix.R
       """
    }
}

(circrna_matrix_mature_seq, circrna_matrix_parent_gene, circrna_matrix_diff_exp) = circRNA_counts.into(3)

process remove_unwanted_biotypes{

    input:
        file(gtf) from ch_gencode_gtf

    output:
        file("filt.gtf") into ch_gtf_filtered

    when: 'circrna_discovery' in module

    script:
    """
    cp ${projectDir}/bin/unwanted_biotypes.txt ./

    grep -vf unwanted_biotypes.txt $gtf > filt.gtf
    """
}

process get_mature_seq{

    publishDir "${params.outdir}/circrna_discovery", mode: params.publish_dir_mode, pattern: 'bed12/*.bed'
    publishDir "${params.outdir}/circrna_discovery", mode: params.publish_dir_mode, pattern: 'fasta/*.fa'

    input:
        file(fasta) from ch_fasta
        file(fai) from ch_fai
        file(gtf) from ch_gtf_filtered
        file(circRNA) from circrna_matrix_mature_seq

    output:
        file("miranda/*.fa") into miranda_sequences
        file("targetscan/*.txt") into targetscan_sequences
        file("bed12/*.bed") into bed_files
        file("fasta/*.fa") into circ_seqs

    when: 'circrna_discovery' in module

    script:
    """
    # convert circrna matrix to bed6 file
    tail -n +2 circRNA_matrix.txt | awk '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, "0", \$4}' | tr ' ' '\t' > circs.bed

    # Create BED12 files
    bash ${projectDir}/bin/get_mature_seq.sh

    # Create miRanda inputs
    bedtools getfasta -fi $fasta -bed de_circ_exon_annotated.bed -s -split -name > de_circ_sequences.fa_tmp
    grep -A 1 '>' de_circ_sequences.fa_tmp | cut -d: -f1,2,3 > de_circ_sequences.fa && rm de_circ_sequences.fa_tmp
    mkdir -p miranda
    awk -F '>' '/^>/ {F=sprintf("miranda/%s.fa",\$2); print > F;next;} {print >> F;}' < de_circ_sequences.fa

    # Create TargetScan inputs
    bedtools getfasta -fi $fasta -bed de_circ_exon_annotated.bed -s -split -tab | sed 's/(/:/g' | sed 's/)//g' > de_circ_seq_tab.txt_tmp
    awk -v OFS="\t" '{print \$1, 9606, \$2}' de_circ_seq_tab.txt_tmp > de_circ_seq_tab.txt && rm de_circ_seq_tab.txt_tmp
    mkdir -p targetscan
    while IFS='' read -r line; do name=\$(echo \$line | awk '{print \$1}'); echo \$line | sed 's/ /\t/g' >> targetscan/\${name}.txt; done < de_circ_seq_tab.txt

    # Save fasta sequences for users
    cp -r miranda/ fasta/
    """
}

(fasta_mature_len, fasta_miranda) = miranda_sequences.into(2)

process get_parent_gene{

    input:
        file(gtf) from ch_gtf_filtered
        file(circRNA) from circrna_matrix_parent_gene

    output:
        file("parent_genes/*.txt") into parent_genes

    when: 'circrna_discovery' in module

    script:
    """
    # convert circrna matrix to bed6 file
    tail -n +2 circRNA_matrix.txt | awk '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, "0", \$4}' | tr ' ' '\t' > circs.bed

    bash ${projectDir}/bin/get_parent_genes.sh
    """
}

process get_mature_len{

    input:
        file(fasta) from fasta_mature_len.flatten()

    output:
        file("*.mature_len.txt") into mature_len

    when: 'circrna_discovery' in module

    script:
    prefix = fasta.toString() - ~/.fa/
    """
    grep -v '>' $fasta | wc -c > ${prefix}.mature_len.txt
    """
}

// Create tuples, merge channels by simpleName for annotation.
ch_mature_len = mature_len.map{ file -> [file.simpleName, file]}
ch_parent_genes = parent_genes.flatten().map{ file -> [file.simpleName, file]}
ch_bed = bed_files.flatten().map{ file -> [file.simpleName, file]}

(bed_ann, bed_circos, bed_diff_exp) = ch_bed.into(3)
(mature_ann, mature_circos, mature_diff_exp) = ch_mature_len.into(3)
(parent_ann, parent_circos, parent_diff_exp) = ch_parent_genes.into(3)

ch_annotate = bed_ann.join(mature_ann).join(parent_ann)

process annotate_circrnas{

    input:
        tuple val(base), file(bed), file(mature_length), file(parent_gene) from ch_annotate

    output:
        file("*annotated.txt") into circrna_annotated

    when: 'circrna_discovery' in module

    script:
    """
    Rscript ${projectDir}/bin/annotate_circs.R $parent_gene $bed $mature_length
    """
}

process master_annotate{

    publishDir "${params.outdir}/circrna_discovery/annotated", mode: params.publish_dir_mode


    input:
        file(annotated) from circrna_annotated.collect()

    output:
        file("circrnas_annotated.txt") into annotated_merged

    when: 'circrna_discovery' in module

    script:
    """
    cat *.txt > merged.txt
    grep -v "Type" merged.txt > no_headers.txt
    echo "circRNA_ID Type Mature_Length Parent_Gene Strand" | tr ' ' '\t' > headers.txt
    cat headers.txt no_headers.txt > circrnas_annotated.txt
    """
}

/*
================================================================================
                         circRNA - miRNA prediction
================================================================================
*/

process miRanda{

    publishDir "${params.outdir}/mirna_prediction/miranda", pattern: "*.miRanda.txt", mode: params.publish_dir_mode

    input:
        file(mirbase) from miranda_miRs
        file(miranda) from fasta_miranda.flatten()

    output:
        file("*.miRanda.txt") into miranda_out

    when: 'mirna_prediction' in module

    script:
    prefix = miranda.toString() - ~/.fa/
    """
    miranda $mirbase $miranda -strict -out ${prefix}.bindsites.out -quiet
    echo "miRNA Target Score Energy_KcalMol Query_Start Query_End Subject_Start Subject_End Aln_len Subject_Identity Query_Identity" | tr ' ' '\t' > ${prefix}.miRanda.txt
    grep -A 1 "Scores for this hit:" ${prefix}.bindsites.out | sort | grep ">" | cut -c 2- | tr ' ' '\t' >> ${prefix}.miRanda.txt
    """
}

process targetscan{

    publishDir "${params.outdir}/mirna_prediction/targetscan", mode: params.publish_dir_mode

    input:
        file(miR) from targetscan_miRs
        file(circ) from targetscan_sequences.flatten()

    output:
        file("*.targetscan.txt") into targetscan_out

    when: 'mirna_prediction' in module

    script:
    prefix = circ.toString() - ~/.txt/
    """
    targetscan_70.pl $miR $circ ${prefix}.targetscan.txt
    """
}

// Create tuples, merge channels by simpleName for report.
ch_targetscan = targetscan_out.map{ file -> [file.simpleName, file]}
ch_miranda = miranda_out.map{ file -> [file.simpleName, file]}

(targetscan_circos, targetscan_diff_exp) = ch_targetscan.into(2)
(miranda_circos, miranda_diff_exp) = ch_miranda.into(2)

ch_circos_plot = targetscan_circos.join(miranda_circos).join(bed_circos).join(parent_circos).join(mature_circos)

process circos_plots{

    publishDir "${params.outdir}/mirna_prediction/circos_plots", pattern: "*.pdf", mode: params.publish_dir_mode
    publishDir "${params.outdir}/mirna_prediction/mirna_targets", pattern: "*miRNA_targets.txt", mode: params.publish_dir_mode

    input:
        tuple val(base), file(targetscan), file(miranda), file(bed), file(parent_gene), file(mature_length) from ch_circos_plot

    output:
        file("*.pdf") into circos_plots
        file("*miRNA_targets.txt") into circrna_mirna_targets

    when: 'mirna_prediction' in module

    script:
    """
    # create file for circos plot
    bash ${projectDir}/bin/prep_circos.sh $bed

    # remove 6mers from TargetScan
    #grep -v "6mer" $targetscan > targetscan_filt.txt

    # Make plots and generate circRNA info
    Rscript ${projectDir}/bin/mirna_circos.R $parent_gene $bed $miranda targetscan_filt.txt $mature_length circlize_exons.txt
    """
}


/*
================================================================================
                          Differential Expression
================================================================================
*/

/*
 * RNA-Seq quantification required to estimate RNA size factors
 * for circRNA library correction
 */

ch_hisat2_index_files = params.hisat2_index ? Channel.value(file(params.hisat2_index + "/*")) : hisat2_built

process Hisat2_align{

    input:
        tuple val(base), file(fastq) from hisat2_reads
        file(hisat2_index) from ch_hisat2_index_files.collect()
        file(fasta) from ch_fasta

    output:
        tuple val(base), file("${base}.bam") into hisat2_bam

    when: 'differential_expression' in module

    script:
    """
    hisat2 -p ${task.cpus} --dta -q -x ${fasta.baseName} -1 ${fastq[0]} -2 ${fastq[1]} -t | samtools view -bS - | samtools sort --threads ${task.cpus} -m 2G - > ${base}.bam
    """
}


process StringTie{

    input:
        tuple val(base), file(bam) from hisat2_bam
        file(gtf) from ch_gencode_gtf

    output:
        file("${base}") into stringtie_dir

    when: 'differential_expression' in module

    script:
    """
    mkdir ${base}/
    stringtie $bam -e -G $gtf -C ${base}/${base}_cov.gtf -p ${task.cpus} -o ${base}/${base}.gtf -A ${base}/${base}_genes.list
    """
}

ch_phenotype = params.phenotype ? file(params.phenotype) : ''

process diff_exp{

    publishDir "${params.outdir}/differential_expression", pattern: "circRNA", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "RNA-Seq", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "boxplots", mode: params.publish_dir_mode
    publishDir "${params.outdir}/quality_control", pattern: "DESeq2_QC", mode: params.publish_dir_mode

    input:
        file(gtf_dir) from stringtie_dir.collect()
        file(circ_matrix) from circrna_matrix_diff_exp
        file(phenotype) from ch_phenotype

    output:
        file("RNA-Seq") into rnaseq_dir
        file("circRNA") into circrna_dir
        file("boxplots") into boxplots_dir
        file("DESeq2_QC") into qc_plots

    when: 'differential_expression' in module

    script:
    """
    for i in \$(ls -d */); do sample=\${i%"/"}; file=\${sample}.gtf; touch samples.txt; printf "\$sample\t\${i}\${file}\n" >> samples.txt; done

    prepDE.py -i samples.txt

    Rscript ${projectDir}/bin/DEA.R gene_count_matrix.csv $phenotype $circ_matrix
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
        'uroborus'
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

// Define STAR index files
def defineStarFiles() {
    return [
    'chrLength.txt',
    'chrNameLength.txt',
    'chrName.txt',
    'chrStart.txt',
    'exonGeTrInfo.tab',
    'exonInfo.tab',
    'geneInfo.tab',
    'Genome',
    'genomeParameters.txt',
    'Log.out',
    'SA',
    'SAindex',
    'sjdbInfo.txt',
    'sjdbList.fromGTF.out.tab',
    'sjdbList.out.tab',
    'transcriptInfo.tab'
    ]
}

// Define Genome versions
def defineGenomeVersions() {
    return [
    'GRCh37',
    'GRCh38'
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
           [ samples, [read1, read2] ]
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
               .filter{ it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
               .ifEmpty{exit 1, "[nf-core/circrna] error: Your FASTQ files do not have the appropriate extension of either '.fastq.gz', 'fq.gz', fastq' or 'fq'."}
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

/*
 * User can supply many explanatory variables and thus this function only checks the condition column.
 * Ask @nf-core how this can be adapted to accept any number of columns?
 */

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

        return condition

        }
        .toList()
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
if(params.skip_trim == 'no'){
    ch_multiqc_report = multiqc_trim_out
}else{
    ch_multiqc_report = multiqc_raw_out
}

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
                log.warn "[nf-core/chipseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/chipseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
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
                    log.error '====================================================\n' +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            '============================================================'
                }
            }
        }
    }
}
