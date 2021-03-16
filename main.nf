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

      --module                        [str] Specify the analysis module(s) to run as a comma seperated list.
                                            Available: 'circrna_discovery', 'mirna_prediction', 'differential_expression'.
                                            Please note that 'circrna_disocvery' is mandatory.

      --tool                          [str] Specify which circRNA quantification tool to use.
                                            If selecting multiple tools, provide as a comma seperated list.
                                            Available: 'circexplorer2', 'circrna_finder', 'ciriquant', 'dcc', 'find_circ', 'mapsplice'.

      --input_type                    [str] Specify the input data type. Avaialable: 'fastq' or 'bam'.

    Input files
      --input                        [file] Either paths to FASTQ/BAM data (must be surrounded with quotes).
                                            Indicate multiple files with a suitable wildcard glob pattern i.e '*_{1,2}' for FASTQ paired end reads.

                                            OR

                                            A path to a CSV file (ending .csv) containing file URL/paths and sample IDs.
                                            Please see documentation for a template.

      --phenotype                    [file] When 'differential_expression' is provided to '--module', a phenotype CSV file (ending .csv) must be provided.
                                            Sample IDs and response + explanatory variables from metadata must be included. Do not include irrelavant metadata.
                                            This file is used to construct the DESeq2 model design.
                                            The response variable must be named 'condition' &  wild-type/control/normal samples named 'control'.
                                            Please see online documentation for examples of a valid phenotype.csv file.

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
 * The below parameters are allowed to be empty (they will be generated if empty)
 * Mainly concerned about valid file extensions when provided.
 */

// Check Fasta
if(params.fasta && (!has_extension(params.fasta, ".fa") || !has_extension(params.fasta, ".fasta"))){
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

/*
 * End of 'allowed to be empty' files
 */

// Check phenotype file

// Check it is a valid file & stage path:
pheno_path = null
if('differential_expression' in module && params.phenotype && (has_extension(params.phenotype, ".csv"))){
	pheno_path = params.phenotype
}else{
	exit 1, "[nf-core/circrna] error: Attempting to run differential expression analysis but the input phenotype file is incorrect.\n\nMust be a '.csv' file and be comma delimited. See online documentation for description + examples."
}

// Check 'condition' is a col name, and that it contains 'control'.
ch_phenotype = Channel.empty()
if(pheno_path){

        pheno_file = file(pheno_path)

	if(pheno_file instanceof List) exit 1, "[nf-core/circrna] error: multiple files passed to --phenotype parameter."
        if(!pheno_file.exists()) exit 1, "[nf-core/circrna] error: input phenotype file could not be found using the path provided: ${params.phenotype}"

        ch_phenotype = examine_phenotype(pheno_file)
	ch_phenotype.filter{ it =~/control/ }
		    .ifEmpty{ exit 1, "[nf-core/circrna] error: There are no samples named control in your condition column. Please rename the wild-type/normal samples control"}

}

// Check BBDUK params

// Check adapters
if(params.skip_trim == 'no'){
	if(!has_extension(params.adapters, ".fa") && !has_extension(params.adapters, ".fasta")){
		exit 1, "[nf-core/circrna] error: --adapters file provied (${params.adapters}) must be a fasta file."
	}
	if(params.adapters){
	adapters = file(params.adapters, checkIfExists: true)
	}
}

// Check all adapter trimming flags are provided
if(params.skip_trim == 'no' && params.adapters && (!params.k && !params.ktrim || !params.k && params.ktrim || params.k && !params.ktrim)){
  exit 1, "[nf-core/circrna] error: Adapter file provided for trimming but missing values for '--k' and/or '--ktrim'.\n\nPlease check the parameter documentation online."
}

// Check all quality trimming flags are provided
if(params.skip_trim == 'no' && (params.trimq && !params.qtrim || !params.trimq && params.qtrim)){
  exit 1, "[nf-core/circrna] error: Both '--trimq' and '--qtrim' are required to perform quality filtering - only one has been provided.\n\nPlease check the parameter documentation online."
}


/*
================================================================================
                          Download Files
================================================================================
*/

process download_genome {

        publishDir "${params.outdir}/circrna_discovery/reference", mode: 'copy'

        output:
          file('*.fa') into fasta_downloaded
          file('*.txt') into gene_annotation_created
          file('*.gtf') into gencode_gtf_downloaded

        when: !(params.fasta) && !(params.gencode_gtf) && !(params.gene_annotation)

        shell:
        if(params.genome_version == 'GRCh37'){
              $/
              wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz
              wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
              gunzip gencode.v34lift37.annotation.gtf.gz
              gunzip GRCh37.primary_assembly.genome.fa.gz
              mv gencode.v34lift37.annotation.gtf GRCh37.gtf
              mv GRCh37.primary_assembly.genome.fa GRCh37.fa.tmp
              sed 's/\s.*$//' GRCh37.fa.tmp > GRCh37.fa
              gtfToGenePred -genePredExt -geneNameAsName2 GRCh37.gtf GRCh37.genepred
              perl -alne '$"="\t";print "@F[11,0..9]"' GRCh37.genepred > GRCh37.txt
              rm GRCh37.fa.tmp
              /$
        }else if(params.genome_version == 'GRCh38'){
              $/
              wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
              wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
              gunzip gencode.v34.primary_assembly.annotation.gtf.gz
              gunzip GRCh38.primary_assembly.genome.fa.gz
              mv gencode.v34.primary_assembly.annotation.gtf GRCh38.gtf
              mv GRCh38.primary_assembly.genome.fa GRCh38.fa.tmp
              sed 's/\s.*$//' GRCh38.fa.tmp > GRCh38.fa
              gtfToGenePred -genePredExt -geneNameAsName2 GRCh38.gtf GRCh38.genepred
              perl -alne '$"="\t";print "@F[11,0..9]"' GRCh38.genepred > GRCh38.txt
              rm GRCh38.fa.tmp
              /$
       }
}

// ternary operators: result = condition ? value_if_true : value_if_false
ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded
ch_gene_annotation = params.gene_annotation ? Channel.value(file(params.gene_annotation)) : gene_annotation_created
ch_gencode_gtf = params.gencode_gtf ? Channel.value(file(params.gencode_gtf)) : gencode_gtf_downloaded

process download_mirbase{
      	errorStrategy 'retry'
        maxRetries 10

      	publishDir "${params.outdir}/mirna_prediction/assets", mode:'copy'

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

// TO DO: add a retry attempt for process below (it sometimes fails to resolve the link)

process download_targetscan{
      	errorStrategy 'retry'
        maxRetries 10

      	publishDir "${params.outdir}/mirna_prediction/assets", mode:'copy'

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
                          Create Genome Index
================================================================================
*/

process samtools_index{

        publishDir "${params.outdir}/circrna_discovery/index/samtools", mode:'copy'

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

        publishDir "${params.outdir}/circrna_discovery/index/bwa", mode:'copy'

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

ch_bwa_index.view()

process hisat2_index{

        publishDir "${params.outdir}/circrna_discovery/index/hisat2", mode: 'copy'

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
ch_hisat2_index.view()

process star_index{

        publishDir "${params.outdir}/circrna_discovery/index", mode:'copy'

        input:
          file(fasta) from ch_fasta
          file(gtf) from ch_gencode_gtf

        output:
          file("STAR") into star_built

        when: !(params.star_index) && ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

        script:
        """
        mkdir star_index

        STAR \
        --runMode genomeGenerate \
        --runThreadN ${params.threads} \
        --sjdbOverhang ${params.sjdbOverhang} \
        --sjdbGTFfile $gtf \
        --genomeDir STAR/ \
        --genomeFastaFiles $fasta
        """
}

ch_star_index = params.star_index ? Channel.value(file(params.star_index)) : star_built
ch_star_index.view()

process bowtie_index{

        publishDir "${params.outdir}/circrna_discovery/index/bowtie", mode:'copy'

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
ch_bowtie_index.view()

process bowtie2_index{

        publishDir "${params.outdir}/circrna_discovery/index/bowtie2", mode:'copy'

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
ch_bowtie2_index.view()


/*
================================================================================
                       Misc. circRNA requirements
================================================================================
*/


process split_fasta{

        publishDir "${params.outdir}/circrna_discovery/reference/chromosomes", mode:'copy'

        input:
          file(fasta) from ch_fasta

        output:
          path("*.fa", includeInputs:true) into split_fasta
          val("${launchDir}/${params.outdir}/circrna_discovery/reference/chromosomes") into split_fasta_path

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

ch_fasta_chr = params.fasta_chr ? Channel.value(params.fasta_chr) : split_fasta_path
ch_fasta_chr.view()

process ciriquant_yml{

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/ciriquant", mode:'copy'

        input:
          file(gencode_gtf) from ch_gencode_gtf
          file(fasta) from ch_fasta
          val(bwa_path) from ch_bwa_index
          val(hisat2_path) from ch_hisat2_index

        output:
          file("travis.yml") into yml_built

        when: !(params.ciriquant_yml) && 'ciriquant' in tool && 'circrna_discovery' in module

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

ch_ciriquant_yml = params.ciriquant_yml ? Channel.value(file(params.ciriquant_yml)) : yml_built

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
if( csv_path ){

    csv_file = file(csv_path)
    println csv_file
    if(csv_file instanceof List) exit 1, "[nf-core/circrna] error: can only accept one csv file per run."
    if(!csv_file.exists()) exit 1, "[nf-core/circrna] error: input CSV file could not be found using the path provided: ${params.input}"

    ch_input_sample = extract_data(csv_path)

}else if (params.input && !has_extension(params.input, "csv")){

    log.info ""
    log.info "No CSV file provided, reading input files from directory provided"
    log.info "Reading path: ${params.input}\n"
    inputSample = retrieve_input_paths(params.input, params.input_type)
    ch_input_sample = inputSample

}else exit 1, "[nf-core/circrna] error: --input file(s) not correctly supplied or improperly defined, see '--help' flag or documentation"

// Process bam file / stage fastq

if(params.input_type == 'bam'){

        process bam_to_fq{

                publishDir "${params.outdir}/quality_control/preprocessing/bamtofastq", mode:'copy'

                input:
                  tuple val(base), file(bam) from ch_input_sample

                output:
                  tuple val(base), file('*.fq.gz') into fastq_built

                script:
                """
                picard -Xmx8g \
                SamToFastq \
                I=$bam \
                F=${base}_R1.fq.gz \
                F2=${base}_R2.fq.gz \
                VALIDATION_STRINGENCY=LENIENT
                """
   }

        (fastqc_reads, trimming_reads, raw_reads, check_reads) = fastq_built.into(4)

}else if(params.input_type == 'fastq'){

        (fastqc_reads, trimming_reads, raw_reads, check_reads) = ch_input_sample.into(4)

}else exit 1, "[nf-core/circrna] error: --input_type must be one of 'fastq' or 'bam'."

check_reads.view()

// FASTQC on raw data. Mandatory.

process FastQC {

        label 'py3'

        publishDir "${params.outdir}/quality_control/fastqc/raw", mode:'copy'

        input:
          tuple val(base), file(fastq) from fastqc_reads

        output:
          file("*.{html,zip}") into fastqc_raw

        script:
        """
        fastqc -q $fastq
        """
}

// MultiQC of the Raw Data, Mandatory.

process multiqc_raw {

	      publishDir "${params.outdir}/quality_control/multiqc", mode:'copy'

      	label 'py3'

      	input:
      	file(htmls) from fastqc_raw.collect()

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

                publishDir "${params.outdir}/quality_control/preprocessing/BBDUK", pattern: "*fq.gz", mode: 'copy'

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
                bbduk.sh -Xmx4g \
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

                publishDir "${params.outdir}/quality_control/fastqc/trimmed", mode:'copy'

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

              	publishDir "${params.outdir}/quality_control/multiqc", mode:'copy'

              	label 'py3'

              	input:
              	  file(htmls) from fastqc_trimmed.collect()
                  file(bbduk_stats) from bbduk_stats_ch.collect()

              	output:
              	  file("Trimmed_Reads_MultiQC.html") into multiqc_trim_out

              	script:
              	"""
              	multiqc -i "Trimmed_Reads_MultiQC" -b "nf-circ pipeline" -n "Trimmed_Reads_MultiQC.html" .
              	"""
 }

}else if(params.skip_trim == 'yes'){
        aligner_reads = raw_reads
}else{
  exit 1, "[nf-core/circrna] error: --skip_trim not specified, please select 'true' or 'false'. See '--help' or documentation for help."
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

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/STAR/1st_Pass", pattern: "${base}", mode:'copy'

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
        --runThreadN ${params.threads} \
        --sjdbScore ${params.sjdbScore} \
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
        """
}

process sjdbFile{

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/STAR/SJFile", pattern: "*SJFile.tab", mode:'copy'

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

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/STAR/2nd_Pass", pattern: "${base}", mode: 'copy'

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
        --runThreadN ${params.threads} \
        --sjdbFileChrStartEnd ${sjdbfile} \
        --sjdbScore ${params.sjdbScore} \
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
        """
}

// CIRCexplorer2

process circexplorer2_star{

        publishDir "${params.outdir}/circrna_discovery/filtered_outputs/circexplorer2", pattern: "*_circexplorer2.bed", mode:'copy'
        publishDir "${params.outdir}/circrna_discovery/tool_outputs/circexplorer2", pattern: "${base}", mode:'copy'

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

        awk '{if(\$13 > 1) print \$0}' ${base}/${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_circexplorer2.bed
        """
}

// circRNA_finder

process circrna_finder{

        publishDir "${params.outdir}/circrna_discovery/filtered_outputs/circrna_finder", pattern: '*_circrna_finder.bed', mode:'copy'
        publishDir "${params.outdir}/circrna_discovery/tool_outputs/circrna_finder/${base}", pattern: "{*filteredJunctions*,*.Chimeric.out.sorted.*}", mode:'copy'

        input:
          tuple val(base), file(star_dir) from circrna_finder_input

        output:
          tuple val(base), file("${base}_circrna_finder.bed") into circrna_finder_results
          tuple val(base), file("{*filteredJunctions*,*.Chimeric.out.sorted.*}") into circrna_finder_raw

        when: 'circrna_finder' in tool && 'circrna_discovery' in module

        script:
        """
        postProcessStarAlignment.pl --starDir ${star_dir}/ --outDir ./

	      awk '{if(\$5 > 1) print \$0}' ${base}.filteredJunctions.bed | awk  -v OFS="\t" -F"\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_circrna_finder.bed
        """
}

// DCC

process dcc_mate1{

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/dcc/${base}", pattern: "mate1", mode:'copy'

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
        --runThreadN ${params.threads} \
        --sjdbFileChrStartEnd ${sjdbfile} \
        --sjdbScore ${params.sjdbScore} \
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
	      """
}

process dcc_mate2{

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/dcc/${base}", pattern: "mate2", mode:'copy'

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
        --runThreadN ${params.threads} \
        --sjdbFileChrStartEnd ${sjdbfile} \
        --sjdbScore ${params.sjdbScore} \
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
	      """
}

ch_dcc_dirs = dcc_pairs.join(dcc_mate1).join(dcc_mate2)

process dcc{

      	label 'py3'

        publishDir "${params.outdir}/circrna_discovery/filtered_outputs/dcc", pattern: "${base}_dcc.bed", mode:'copy'
        publishDir "${params.outdir}/circrna_discovery/tool_outputs/dcc/${base}", pattern: "{*.log,*Circ*}", mode:'copy'

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
        DCC @samplesheet -mt1 @mate1file -mt2 @mate2file -D -an $gtf -Pi -ss -F -M -Nr 1 1 -fg -A $fasta -N -T ${params.threads}
        awk '{print \$6}' CircCoordinates >> strand
        paste CircRNACount strand | tail -n +2 | awk -v OFS="\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${base}_dcc.txt
	      bash ${projectDir}/bin/filter_DCC.sh ${base}_dcc.txt

        mv CircCoordinates ${base}.CircCoordinates
        mv CircRNACount ${base}.CircRNACount
        mv *.log ${base}.log
        """
}

// find_circ

process find_anchors{

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/find_circ/${base}", pattern: "{*anchors.qfa.gz,*.bam}", mode:'copy'

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
        bowtie2 -p ${params.threads} --very-sensitive --mm -D 20 --score-min=C,-15,0 \
        -x ${fasta.baseName} -q -1 ${fastq[0]} -2 ${fastq[1]} \
        | samtools view -hbuS - | samtools sort --threads ${params.threads} -m 2G - > ${base}.bam

        samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

        unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz
        """
}

process find_circ{

        publishDir "${params.outdir}/circrna_discovery/filtered_outputs/find_circ/", pattern: '*_find_circ.bed', mode:'copy'
        publishDir "${params.outdir}/circrna_discovery/tool_outputs/find_circ/${base}", pattern: "*.sites.*", mode: 'copy'

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
        bowtie2 -p ${params.threads} --reorder --mm -D 20 --score-min=C,-15,0 -q -x ${fasta.baseName} \
        -U $anchors | python ${projectDir}/bin/find_circ.py -G $fasta_chr_path -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

        echo "# chrom:start:end:name:n_reads:strand:n_uniq:best_qual_A:best_qual_B:spliced_at_begin:spliced_at_end:tissues:tiss_counts:edits:anchor_overlap:breakpoints" > tmp.txt

        cat tmp.txt | tr ':' '\t' > ${base}.bed

        grep circ ${base}.sites.bed | grep -v chrM | python ${projectDir}/bin/sum.py -2,3 | python ${projectDir}/bin/scorethresh.py -16 1 | python ${projectDir}/bin/scorethresh.py -15 2 | python ${projectDir}/bin/scorethresh.py -14 2 | python ${projectDir}/bin/scorethresh.py 7 2 | python ${projectDir}/bin/scorethresh.py 8,9 35 | python ${projectDir}/bin/scorethresh.py -17 100000 >> ${base}.txt

	      tail -n +2 ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_find_circ.bed
	      """
}

// CIRIquant

process ciriquant{

        publishDir "${params.outdir}/circrna_discovery/filtered_outputs/ciriquant", pattern: "${base}_ciriquant.bed", mode:'copy'
        publishDir "${params.outdir}/circrna_discovery/tool_outputs/ciriquant", pattern: "${base}", mode: 'copy'

        input:
          tuple val(base), file(fastq) from ciriquant_reads
          file(ciriquant_yml) from ch_ciriquant_yml

        output:
          tuple val(base), file("${base}_ciriquant.bed") into ciriquant_results
          tuple val(base), file("${base}") into ciriquant_raw_dir

        when: 'ciriquant' in tool && 'circrna_discovery' in module

        script:
        """
        CIRIquant -t ${params.threads} \
        -1 ${fastq[0]} \
        -2 ${fastq[1]} \
        --config $ciriquant_yml \
        --no-gene \
        -o ${base} \
        -p ${base}

        cp ${base}/${base}.gtf ${base}_ciriquant.gtf

	      bash filter_CIRIquant.sh ${base}_ciriquant.gtf
        mv ${base}_ciriquant.gtf ${base}.gtf
        """
}


// mapsplice

process mapsplice_align{

        publishDir "${params.outdir}/circrna_discovery/tool_outputs/mapsplice", pattern: "${base}", mode:'copy'

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
            -p ${params.threads} \
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
            -p ${params.threads} \
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

        publishDir "${params.outdir}/circrna_discovery/filtered_outputs/mapsplice", pattern: "*_mapsplice.bed", mode:'copy'

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

	      awk '{if(\$13 > 1) print \$0}' ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_mapsplice.bed
        """
}

// UROBORUS

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
        tophat -p ${params.threads} -o ${base} ${fasta.baseName} ${fastq[0]} ${fastq[1]}
        mv ${base}/unmapped.bam ./
        mv ${base}/accepted_hits.bam ./
        """
}


process uroborus{

        publishDir "${params.outdir}/circrna_discovery/uroborus", mode:'copy'

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

/*
 * CONSOLIDATION OF TOOLS
 * Keep circRNAs that have been called by at least 2 tools (if tools_selected > 1)
 * circRNA tool outputs converted to count matrix to facilitate filtering
 * (circRNAs with read counts < 1 have already been removed during aligner processes)
 */

// check the length of the tool list
tools_selected = tool.size()

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

          ## Bring forward circRNAs called by at least 2 tools
          Rscript ${projectDir}/bin/consolidate_algorithms.R samples.csv

          mv combined_counts.bed ${base}.bed
          """
  }

  process get_counts_combined{

          publishDir "${params.outdir}/circrna_discovery/count_matrix", mode:'copy'

			    input:
				    file(bed) from sample_counts.collect()

			    output:
				    file("circRNA_matrix.txt") into circRNA_counts

          when: 'circrna_discovery' in module

          script:
      		"""
      		python ${projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
      		"""
		}

} else{

  single_tool = ciriquant_results.mix(circexplorer2_results, dcc_results, circrna_finder_results, find_circ_results, mapsplice_results)

  process get_counts_single{

          echo true
          publishDir "${params.outdir}/circrna_discovery/count_matrix", mode:'copy'

          input:
            file(bed) from single_tool.collect()
				    val(tool) from params.tool

          output:
            file("circRNA_matrix.txt") into circRNA_counts

          when: 'circrna_discovery' in module

          script:
          """
          for b in *.bed; do
				      foo=\${b%".bed"};
				      bar=\${foo%"_${tool}"};
				      mv \$b \${bar}.bed
			    done

			    python ${projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
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

        publishDir "${params.outdir}/circrna_discovery", mode:'copy', pattern: 'bed12/*.bed'
        publishDir "${params.outdir}/circrna_discovery", mode:'copy', pattern: 'fasta/*.fa'

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

        publishDir "${params.outdir}/circrna_discovery/annotated", mode: 'copy'

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

	      publishDir "${params.outdir}/mirna_prediction/miranda", pattern: "*.miRanda.txt", mode:'copy'

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

    	  publishDir "${params.outdir}/mirna_prediction/targetscan", mode:'copy'

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

        publishDir "${params.outdir}/mirna_prediction/circos_plots", pattern: "*.pdf", mode: 'copy'
        publishDir "${params.outdir}/mirna_prediction/mirna_targets", pattern: "*miRNA_targets.txt", mode: 'copy'

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
      	grep -v "6mer" $targetscan > targetscan_filt.txt

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
        hisat2 -p 2 --dta -q -x ${fasta.baseName} -1 ${fastq[0]} -2 ${fastq[1]} -t | samtools view -bS - | samtools sort --threads 2 -m 2G - > ${base}.bam
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
        stringtie $bam -e -G $gtf -C ${base}/${base}_cov.gtf -p 2 -o ${base}/${base}.gtf -A ${base}/${base}_genes.list
        """
}

//if('differential_expression' in module && (!params.phenotype || params.phenotype == null)){
//  exit 1, "[nf-core/circrna] error: parameter '--phenotype' (file for DESeq2) not supplied. Please see '--help' or documentation for help."
//}else if('differential_expression' in module && (params.phenotype != null)){
//  ch_phenotype = file(params.phenotype)
//}

ch_phenotype = params.phenotype ? file(params.phenotype) : ''

process diff_exp{

        publishDir "${params.outdir}/differential_expression", pattern: "circRNA", mode:'copy'
        publishDir "${params.outdir}/differential_expression", pattern: "RNA-Seq", mode:'copy'
        publishDir "${params.outdir}/differential_expression", pattern: "boxplots", mode:'copy'
        publishDir "${params.outdir}/quality_control", pattern: "DESeq2_QC", mode:'copy'

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
                         Functions
================================================================================
*/

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
        if (!row.keySet().containsAll(expected_keys)) exit 1, "[nf-core/circrna] error: Invalid CSV input - malformed column names. Please use the column names 'Sample_ID', 'Read1', 'Read2', 'Bam'."

        checkNumberOfItem(row, 4)

        def samples = row.Sample_ID
        def read1 = row.Read1.matches('NA') ? 'NA' : return_file(row.Read1)
        def read2 = row.Read2.matches('NA') ? 'NA' : return_file(row.Read2)
        def bam = row.Bam.matches('NA') ? 'NA' : return_file(row.Bam)

        if(samples == '' || read1 == '' || read2 == '' || bam == '') exit 1, "[nf-core/circrna] error: a field does not contain any information. Please check your CSV file"

        if(read1.matches('NA') && read2.matches('NA') && bam.matches('NA')) exit 1, "[nf-core/circrna] error: A row in your CSV file appears to have missing information."

        if ( !read1.matches('NA') && !has_extension(read1, "fastq.gz") && !has_extension(read1, "fq.gz") && !has_extension(read1, "fastq") && !has_extension(read1, "fq")) exit 1, "[nf-core/circrna] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
        if ( !read2.matches('NA') && !has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz") && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "[nf-core/circrna] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
        if ( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified BAM file either has a non-recognizable extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

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
