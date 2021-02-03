#!/usr/bin/env nextflow

/*
================================================================================
                                circRNA analysis
================================================================================
Started August 2020
--------------------------------------------------------------------------------
Description:
  (To my knowledge) the first circRNA pipeline to scan RNA-Seq data for circRNAs and
  conduct differential expression analysis + circRNA-miRNA predictions
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/BarryDigby/circRNA
 -------------------------------------------------------------------------------
 @Documentation
 Work in progress
--------------------------------------------------------------------------------
*/

/*
================================================================================
                                  Help Flags
================================================================================
*/

ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";


def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }

//help information
params.help = null
if (params.help) {
    log.info ''
    log.info print_purple('------------------------------------------------------------------------')
    log.info "tool_name: A Nextflow based circular RNA analsis pipeline"
    log.info "tool_name integrates several NGS processing tools to identify novel circRNAs from "
    log.info "un-processed RNA sequencing data. To run this pipeline users need to install nextflow"
    log.info "and singularity. "
    print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ') +

            print_purple('Nextflow run BarryDigby/circRNA --profile singularity, standard <options> \n') +

            print_yellow('      Mandatory arguments:\n') +
            print_cyan('      --input <path>            ') + print_green('Path to input data\n') +
            print_cyan('      --input_type <str>           ') + print_green('Input data type. Supported: fastq, bam\n') +
            print_cyan('      --input_glob <str>           ') + print_green('Glob pattern of input files e.g: \'_R{1,2}.fastq.gz\'\n') +
            print_cyan('      --tool <str>                 ') + print_green('circRNA tool to use for analysis. \n') +
            print_green('                                   Supported: CIRCexplorer2, CIRIquant, find_circ\n') +
            print_green('                                   mapsplice, DCC, circRNA_finder\n') +
            print_cyan('      --genome_version <str>              ') + print_green('Supported: GRCh37, GRCh38\n') +
            '\n' +
            print_yellow('    Input Files:            if left empty these will be generated\n') +
            print_cyan('      --fasta <path>               ') + print_green('Path to genome fasta file\n') +
            print_cyan('      --fasta_fai <path>           ') + print_green('Path to genome fasta fai file\n') +
            print_cyan('      --gencode_gtf <path>         ') + print_green('Path to genocde gtf file\n') +
            print_cyan('      --gene_annotation <path>     ') + print_green('Path to gene annotation file \n') +
            print_cyan('      --star_index <str>           ') + print_green('Path to STAR index\n') +
            print_cyan('      --bwa_index <str>            ') + print_green('Path to BWA index\n') +
            print_cyan('      --bowtie_index <str>         ') + print_green('Path to Bowtie index (must include glob for files)\n') +
            print_cyan('      --bowtie2_index <str>        ') + print_green('Path to Bowtie2 index (must include glob for files)\n') +
            print_cyan('      --hisat2_index <str>         ') + print_green('Path to Hisat2 index\n') +
            print_cyan('      --ciriquant_yml <str>        ') + print_green('Path to CIRIquant yml configuration file\n') +
            print_cyan('      --adapters <path>            ') + print_green('Fasta file containing adapters to trim\n')


            log.info ('------------------------------------------------------------------------')
            log.info print_yellow('Contact information: b.digby237@gmail.com')
            log.info print_yellow('O\'Broin Lab, National University of Ireland Galway')
    log.info ('------------------------------------------------------------------------')
    exit 0
}


/*
================================================================================
                                  Paramaters
================================================================================
*/


params.outdir = null
params.fasta = null
params.gencode_gtf = null
params.gene_annotation = null
params.genome_version = null
params.tool = null
params.fasta_fai = null
params.bwa_index = null
params.star_index = null
params.hisat2_index = null
params.bowtie_index = null
params.bowtie2_index = null
params.fasta_chr = null
params.ciriquant_yml = null
params.input = null
params.input_type = null
params.input_glob = null
params.adapters = null
params.phenotype = null
params.trimming = null
params.k = null
params.ktrim = null
params.qtrim = null
params.trimq = null
params.minlen = null

/*
================================================================================
                          Check mandatory flags
================================================================================
*/

// Check Tools selected
toolList = defineToolList()
tool = params.tool ? params.tool.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tool, toolList)) exit 1, "[nf-core/circrna] error: Unknown tool, see --help for more information."

// Check Modules
moduleList = defineModuleList()
module = params.module ? params.module.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(module, moduleList)) exit 1, "[nf-core/circrna] error: Unknown module selected, see --help for more information."

// Check Input parameter
if(params.input == null){
  exit 1, "[nf-core/circrna] error: --input was not supplied! Please check '--help' or documentation for details"
}

// Check input type
if(params.input_type == null){
  exit 1, "[nf-core/circrna] error: --input_type was not supplied! Please select 'fastq' or 'bam'."
}

// Check Genome version
if(params.genome_version == null){
  exit 1, "[nf-core/circrna] error: --genome_version was not supplied!. Please select 'GRCh37' or 'GRCh37'."
}


// Check

/*
================================================================================
                          Download Files
================================================================================
*/

process download_genome {

        publishDir "$params.outdir/reference", mode: 'copy'

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
              mv gencode.v34.primary_assembly.annotation.gtf GRCh37.gtf
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

      	publishDir "$params.outdir/assets", mode:'copy'

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

      	publishDir "$params.outdir/assets", mode:'copy'

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

        publishDir "$params.outdir/reference", mode:'copy'

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

        publishDir "$params.outdir/index/bwa", mode:'copy'

        input:
          file(fasta) from ch_fasta

        output:
          file("${fasta.baseName}.*") into bwa_built
          val("$launchDir/index/bwa") into bwa_path

        when: !(params.bwa_index) && 'ciriquant' in tool && 'circrna_discovery' in module

        script:
        """
        bwa index -a bwtsw $fasta -p ${fasta.baseName}
        """
}

ch_bwa_index = params.bwa_index ? Channel.value(params.bwa_index) : bwa_path

ch_bwa_index.view()

process hisat2_index{

        publishDir "$params.outdir/index/hisat2", mode: 'copy'

        input:
          file(fasta) from ch_fasta

        output:
          file("${fasta.baseName}.*.ht2") into hisat2_built
          val("$launchDir/index/hisat2") into hisat2_path

        when: !(params.hisat2_index) && 'ciriquant' in tool && ('circrna_discovery' || 'differential_expression' in module)

        script:
        """
        hisat2-build $fasta ${fasta.baseName}
        """
}

ch_hisat2_index = params.hisat2_index ? Channel.value(params.hisat2_index) : hisat2_path
ch_hisat2_index.view()

process star_index{

        publishDir "$params.outdir/index", mode:'copy'

        input:
          file(fasta) from ch_fasta
          file(gtf) from ch_gencode_gtf

        output:
          file("star_index") into star_built

        when: !(params.star_index) && ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool) && 'circrna_discovery' in module

        script:
        """
        mkdir star_index

        STAR \
        --runMode genomeGenerate \
        --runThreadN 8 \
        --sjdbGTFfile $gtf \
        --genomeDir star_index/ \
        --genomeFastaFiles $fasta
        """
}

ch_star_index = params.star_index ? Channel.value(file(params.star_index)) : star_built
ch_star_index.view()

process bowtie_index{

        publishDir "$params.outdir/index/bowtie", mode:'copy'

        input:
          file(fasta) from ch_fasta

        output:
          file ("${fasta.baseName}.*") into bowtie_built
          val("$launchDir/index/bowtie") into bowtie_path

        when: !(params.bowtie_index) && ('mapsplice' in tool || 'uroborus' in tool) && 'circrna_discovery' in module

        script:
        """
        bowtie-build $fasta ${fasta.baseName}
        """
}

ch_bowtie_index = params.bowtie_index ? Channel.value(file(params.bowtie_index)) : bowtie_built
ch_bowtie_index.view()

process bowtie2_index{

        publishDir "$params.outdir/index/bowtie2", mode:'copy'

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

ch_bowtie2_index = params.bowtie2_index ? Channel.value(file(params.bowtie2_index)) : bowtie2_built
ch_bowtie2_index.view()


/*
================================================================================
                       Misc. circRNA requirements
================================================================================
*/


process split_fasta{

        publishDir "$params.outdir/index/chromosomes", mode:'copy'

        input:
          file(fasta) from ch_fasta

        output:
          path("*.fa", includeInputs:true) into split_fasta
          val("$launchDir/index/chromosomes") into split_fasta_path

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

        publishDir "$params.outdir", mode:'copy'

        input:
          file(fasta) from ch_fasta
          val(gencode_gtf_path) from ch_gencode_gtf
          val(fasta_path) from ch_fasta
          val(bwa_path) from ch_bwa_index
          val(hisat2_path) from ch_hisat2_index

        output:
          file("travis.yml") into yml_built

        when: !(params.ciriquant_yml) && 'ciriquant' in tool && 'circrna_discovery' in tool

        script:
        index_prefix = fasta.toString() - ~/.fa/
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

    if(csv_file instanceof List) exit 1, "[nf-core/circrna] error: can only accept one csv file per run."
    if(!csv_file.exists()) exit 1, "[nf-core/circrna] error: input CSV file could not be found using the path provided: ${params.input}"

    ch_input_sample = extract_data(csv_path)

}else if (params.input && !has_extension(params.input, "csv")){

    log.info ""
    log.info "No CSV file provided, reading input files from directory provided"
    log.info "Reading path: ${params.input}\n"
    inputSample = retrieve_input_paths(params.input, params.input_glob, params.input_type)
    ch_input_sample = inputSample

}else exit 1, "[nf-core/circrna] error: --input file(s) not correctly supplied or improperly defined, see '--help' flag or documentation"

// Process bam file / stage fastq

if(params.input_type == 'bam'){

        process bam_to_fq{

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

// FASTQC on raw data. Mandatory.

process FastQC {

        publishDir "$params.outdir/quality_control/fastqc/raw", mode:'copy'

        input:
          tuple val(base), file(fastq) from fastqc_reads

        output:
          file("*.{html,zip}") into fastqc_raw

        script:
        """
        fastqc -q $fastq
        """
}

// BBDUK

if(params.trimming == true){

        process bbduk {

                publishDir "$params.outdir/trimmed_reads", mode:'copy'

                input:
                  tuple val(base), file(fastq) from trimming_reads
                  path adapters from params.adapters

                output:
                  tuple val(base), file('*.trim.fq.gz') into trim_reads_ch

                script:
                """
                bbduk.sh -Xmx4g \
                in1=${fastq[0]} \
                in2=${fastq[1]} \
                out1=${base}_1.trim.fq.gz \
                out2=${base}_2.trim.fq.gz \
                ref=$adapters \
                minlen=30 \
                ktrim=r \
                k=12 \
                qtrim=r \
                trimq=20
                """
        }

        // trimmed reads into 2 channels:
        (fastqc_trim_reads, aligner_reads) = trim_reads_ch.into(2)

        process FastQC_trim {

                publishDir "$params.outdir/quality_control/fastqc/trimmed", mode:'copy'

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

              	publishDir "$params.outdir/quality_control/multiqc/trimmed", mode:'copy'

              	label 'multiqc'

              	input:
              	  file(htmls) from fastqc_trimmed.collect()

              	output:
              	  file("Trimmed_Reads_MultiQC.html") into multiqc_trim_out

              	script:
              	"""
              	multiqc -i "Trimmed_Reads_MultiQC" -b "nf-circ pipeline" -n "Trimmed_Reads_MultiQC.html" .
              	"""
 }

}else if(params.trimming == false){
        aligner_reads = raw_reads
}

// MultiQC of the Raw Data, Mandatory.

process multiqc_raw {

	      publishDir "$params.outdir/quality_control/multiqc/raw", mode:'copy'

      	label 'multiqc'

      	input:
      	file(htmls) from fastqc_raw.collect()

      	output:
      	file("Raw_Reads_MultiQC.html") into multiqc_raw_out

      	script:
      	"""
      	multiqc -i "Raw_Reads_MultiQC" -b "nf-circ pipeline" -n "Raw_Reads_MultiQC.html" .
      	"""
}


// Stage Aligner read channels
(circexplorer2_reads, find_circ_reads, ciriquant_reads, mapsplice_reads, uroborus_reads, circrna_finder_reads, dcc_reads, dcc_reads_mate1, dcc_reads_mate2, hisat2_reads) = aligner_reads.into(10)

/*
================================================================================
                             circRNA Discovery
================================================================================
*/

// CIRCexplorer2

process star_align{

        publishDir "$params.outdir/star_alignment", mode:'copy', overwrite: true

        input:
          tuple val(base), file(fastq) from circexplorer2_reads
          file(gtf) from ch_gencode_gtf
          val(star_idx) from ch_star_index

        output:
          tuple val(base), file("${base}.Chimeric.out.junction") into circexplorer2_input

        when: 'circexplorer2' in tool && 'circrna_discovery' in module

        script:
        """
        STAR    \
        --runThreadN 16 \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --genomeLoad NoSharedMemory \
        --genomeDir $star_idx \
        --readFilesIn ${fastq[0]},${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${base}. \
        --outSJfilterOverhangMin 15 15 15 15 \
        --outFilterMultimapNmax 1 \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMin 1 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 10 \
        --outFilterMismatchNoverLmax 0.05 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 1 \
        --alignSJDBoverhangMin 1 \
        --alignSoftClipAtReferenceEnds No \
        --chimSegmentMin 10 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15 \
        --sjdbGTFfile $gtf  \
        --sjdbScore 2 \
        --chimOutType Junctions \
        --outSAMtype BAM SortedByCoordinate
        """
}


process circexplorer2_star{

        publishDir "$params.outdir/circrna_discovery/circexplorer2/parsed", pattern: '*_circexplorer2.bed', mode:'copy'
        publishDir "$params.outdir/circrna_discovery/circexplorer2/raw", pattern: "${base}.txt", mode:'copy'

        input:
          tuple val(base), file(chimeric_reads) from circexplorer2_input
          file(fasta) from ch_fasta
          file(gene_annotation) from ch_gene_annotation

        output:
          tuple val(base), file("${base}_circexplorer2.bed") into circexplorer2_results
          tuple val(base), file("${base}.txt") into circexplorer2_raw_results

        when: 'circexplorer2' in tool && 'circrna_discovery' in module

        script:
        """
        CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}.STAR.junction.bed
        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.STAR.junction.bed -o ${base}.txt

        awk '{if(\$13 > 1) print \$0}' ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_circexplorer2.bed
        """
}


// find_circ

process find_anchors{

        input:
          tuple val(base), file(fastq) from find_circ_reads
          file(fasta) from ch_fasta
          file(bowtie2_index) from ch_bowtie2_index.collect()

        output:
          tuple val(base), file("${base}_anchors.qfa.gz") into ch_anchors

        when: 'find_circ' in tool && 'circrna_discovery' in module

        script:
        """
        bowtie2 -p 16 --very-sensitive --mm -D 20 --score-min=C,-15,0 \
        -x ${fasta.baseName} -q -1 ${fastq[0]} -2 ${fastq[1]} \
        | samtools view -hbuS - | samtools sort --threads 16 -m 2G - > ${base}.bam

        samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

        unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz
        """
}


process find_circ{

        publishDir "$params.outdir/circrna_discovery/find_circ/parsed", pattern: '*_find_circ.bed', mode:'copy'
        publishDir "$params.outdir/circrna_discovery/find_circ/raw", pattern: "${base}.txt", mode: 'copy'

        input:
          tuple val(base), file(anchors) from ch_anchors
          file(bowtie2_index) from ch_bowtie2_index.collect()
          file(fasta) from ch_fasta
          val(fasta_chr_path) from ch_fasta_chr

        output:
          tuple val(base), file("${base}_find_circ.bed") into find_circ_results
          tuple val(base), file("${base}.txt") into find_circ_raw_results

        when: 'find_circ' in tool && 'circrna_discovery' in module

        script:
        """
        bowtie2 -p 16 --reorder --mm -D 20 --score-min=C,-15,0 -q -x ${fasta.baseName} \
        -U $anchors | python /opt/conda/envs/circrna/bin/find_circ.py -G $fasta_chr_path -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

        echo "# chrom:start:end:name:n_reads:strand:n_uniq:best_qual_A:best_qual_B:spliced_at_begin:spliced_at_end:tissues:tiss_counts:edits:anchor_overlap:breakpoints" > tmp.txt

        cat tmp.txt | tr ':' '\t' > ${base}.bed

        grep circ ${base}.sites.bed | grep -v chrM | python /opt/conda/envs/circrna/bin/sum.py -2,3 | python /opt/conda/envs/circrna/bin/scorethresh.py -16 1 | python /opt/conda/envs/circrna/bin/scorethresh.py -15 2 | python /opt/conda/envs/circrna/bin/scorethresh.py -14 2 | python /opt/conda/envs/circrna/bin/scorethresh.py 7 2 | python /opt/conda/envs/circrna/bin/scorethresh.py 8,9 35 | python /opt/conda/envs/circrna/bin/scorethresh.py -17 100000 >> ${base}.txt

	      tail -n +2 ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_find_circ.bed
	      """
}


// circRNA_finder

process circrna_finder_star{

        input:
          tuple val(base), file(fastq) from circrna_finder_reads
          val(star_index) from ch_star_index

        output:
          tuple val(base), file("${base}") into circrna_finder_star

        when: 'circrna_finder' && 'circrna_discovery' in module

        script:
        """
        STAR \
        --genomeDir $star_index \
        --readFilesIn ${fastq[0]} ${fastq[1]} \
        --readFilesCommand zcat \
        --runThreadN 16 \
        --chimSegmentMin 20 \
        --chimScoreMin 1 \
        --chimOutType Junctions SeparateSAMold \
        --alignIntronMax 100000 \
        --outFilterMismatchNmax 4 \
        --alignTranscriptsPerReadNmax 100000 \
        --outFilterMultimapNmax 2 \
        --outFileNamePrefix ${base}/${base}.
        """
}


process circrna_finder{

        publishDir "$params.outdir/circrna_discovery/circrna_finder/parsed", pattern: '*_circrna_finder.bed', mode:'copy'
        publishDir "$params.outdir/circrna_discovery/circrna_finder/raw", pattern: "${base}.filteredJunctions.bed", mode:'copy'

        input:
          tuple val(base), file(star_dir) from circrna_finder_star

        output:
          tuple val(base), file("${base}_circrna_finder.bed") into circrna_finder_results
          tuple val(base), file("${base}.filteredJunctions.bed") into circrna_finder_raw_results

        when: 'circrna_finder' in tool && 'circrna_discovery' in module

        script:
        """
        postProcessStarAlignment.pl --starDir ${star_dir}/ --outDir ./

	      tail -n +2 ${base}.filteredJunctions.bed | awk '{if(\$5 > 1) print \$0}' | awk  -v OFS="\t" -F"\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_circrna_finder.bed
        """
}

// DCC

process dcc_pair{

        input:
          tuple val(base), file(fastq) from dcc_reads
          val(star_index) from ch_star_index

        output:
          tuple val(base), file("samples") into dcc_samples

        when: 'dcc' in tool && 'circrna_discovery' in module

        script:
        """
        STAR \
        --runThreadN 16 \
        --genomeDir $star_index \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${fastq[0]} ${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix samples/${base}. \
        --outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15
        """
}

process dcc_1{

        input:
          tuple val(base), file(fastq) from dcc_reads_mate1
          val(star_index) from ch_star_index

        output:
          tuple val(base), file("mate1") into dcc_mate1

       when: 'dcc' in tool && 'circrna_discovery' in module

        script:
        """
        STAR \
        --runThreadN 16 \
        --genomeDir $star_index \
        --outSAMtype None \
        --readFilesIn ${fastq[0]} \
        --readFilesCommand zcat \
        --outFileNamePrefix mate1/${base}. \
        --outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --seedSearchStartLmax 30 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15
        """
}

process dcc_2{

        input:
          tuple val(base), file(fastq) from dcc_reads_mate2
          val(star_index) from ch_star_index

        output:
          tuple val(base), file("mate2") into dcc_mate2

        when: 'dcc' in tool && 'circrna_discovery' in module

        script:
        """
        STAR \
        --runThreadN 16 \
        --genomeDir $star_index \
        --outSAMtype None \
        --readFilesIn ${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix mate2/${base}. \
        --outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --seedSearchStartLmax 30 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15
        """
}

// collect runs according to val(base) in tuple
ch_dcc_dirs = dcc_samples.join(dcc_mate1).join(dcc_mate2)

process dcc{

        publishDir "$params.outdir/circrna_discovery/dcc/parsed", pattern: "${base}_dcc.txt", mode:'copy'
        publishDir "$params.outdir/circrna_discovery/dcc/raw", pattern: "${base}.Circ*", mode:'copy'

        input:
          tuple val(base), file(samples), file(mate1), file(mate2) from ch_dcc_dirs
          file(gtf) from ch_gencode_gtf
          file(fasta) from ch_fasta

        output:
          tuple val(base), file("${base}_dcc.bed") into dcc_results
          tuple val(base), file("${base}.Circ*") into dcc_raw_results

        when: 'dcc' in tool && 'circrna_discovery' in module

        script:
        COJ="Chimeric.out.junction"
        """
        sed -i 's/^chr//g' $gtf

        printf "samples/${base}.${COJ}" > samplesheet
        printf "mate1/${base}.${COJ}" > mate1file
        printf "mate2/${base}.${COJ}" > mate2file

        DCC @samplesheet -mt1 @mate1file -mt2 @mate2file -D -an $gtf -Pi -ss -F -M -Nr 1 1 -fg -A $fasta -N -T 8

        awk '{print \$6}' CircCoordinates >> strand
        paste CircRNACount strand | tail -n +2 | awk -v OFS="\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${base}_dcc.txt
	      bash filter_DCC.sh ${base}_dcc.txt

        mv CircCoordinates ${base}.CircCoordinates
        mv CircRNACount ${base}.CircRNACount
        """
}

// CIRIquant


process ciriquant{

        publishDir "$params.outdir/circrna_discovery/ciriquant/parsed", pattern: "${base}_ciriquant.bed", mode:'copy'
        publishDir "$params.outdir/circrna_discovery/ciriquant/raw", pattern: "${base}.gtf", mode: 'copy'

        input:
          tuple val(base), file(fastq) from ciriquant_reads
          file(ciriquant_yml) from ch_ciriquant_yml

        output:
          tuple val(base), file("${base}_ciriquant.bed") into ciriquant_results
          tuple val(base), file("${base}.gtf") into ciriquant_raw_results

        when: 'ciriquant' in tool && 'circrna_discovery' in module

        script:
        """
        CIRIquant -t 2 \
        -1 ${fastq[0]} \
        -2 ${fastq[1]} \
        --config $ciriquant_yml \
        --no-gene \
        -o ${base} \
        -p ${base}

        mv ${base}/${base}.gtf ${base}_ciriquant.gtf

	      bash filter_CIRIquant.sh ${base}_ciriquant.gtf
        mv ${base}_ciriquant.gtf ${base}.gtf
        """
}


// mapsplice

process mapsplice_align{

        publishDir "$params.outdir/mapsplice", mode:'copy'

        input:
          tuple val(base), file(fastq) from mapsplice_reads
          val(mapsplice_ref) from ch_fasta_chr
          file(bowtie_index) from ch_bowtie_index.collect()
          file(gtf) from ch_gencode_gtf

        output:
          tuple val(base), file("${base}/fusions_raw.txt") into mapsplice_fusion

        when: 'mapsplice' in tool && 'circrna_discovery' in module

        script:
        prefix = gtf.toString() - ~/.gtf/
        """
        gzip -d --force ${fastq[0]}
        gzip -d --force ${fastq[1]}

        mapsplice.py \
        -c $mapsplice_ref \
        -x $prefix \
        -1 ${base}_r1.fastq \
        -2 ${base}_r2.fastq \
        -p 8 \
        --bam \
        --seglen 25 \
        --min-map-len 40 \
        --fusion-non-canonical \
        --min-fusion-distance 200 \
        --gene-gtf $gtf \
        -o $base
        """
}


process mapsplice_parse{

        publishDir "$params.outdir/circrna_discovery/mapsplice/parsed", pattern: "*_mapsplice.bed", mode:'copy'
        publishDir "$params.outdir/circrna_discovery/mapsplice/raw", pattern: "${base}.txt", mode:'copy'

        input:
          tuple val(base), file(raw_fusion) from mapsplice_fusion
          file(fasta) from ch_fasta
          file(gene_annotation) from ch_gene_annotation

        output:
          tuple val(base), file("${base}_mapsplice.bed") into mapsplice_results
          tuple val(base), file("${base}.txt") into mapsplice_raw_results

        when: 'mapsplice' in tool && 'circrna_discovery' in module

        script:
        """
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
        tophat -p 8 -o ${base} ${fasta.baseName} ${fastq[0]} ${fastq[1]}
        mv ${base}/unmapped.bam ./
        mv ${base}/accepted_hits.bam ./
        """
}


process uroborus{

        publishDir "$params.outdir/circrna_discovery/uroborus", mode:'copy'

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

	combined_tool = ciriquant_results.join(circexplorer2_results).join(dcc_results).join(circrna_finder_results).join(find_circ_results).join(mapsplice_results)

  process consolidate_algorithms{

          echo true
          publishDir "$params.outdir/circrna_discovery/count_matrix", mode:'copy'

          input:
            tuple val(base), file(ciriquant), file(circexplorer2), file(dcc), file(circrna_finder), file(find_circ), file(mapsplice) from combined_tool

          output:
            file("${base}.bed") into sample_counts

          when: ( 'circrna_discovery' || 'differential_expression' || 'mirna_prediction' in module)

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

          publishDir "$params.outdir/circrna_discovery/count_matrix", mode:'copy'

			    input:
				    file(bed) from sample_counts.collect()

			    output:
				    file("circRNA_matrix.txt") into circRNA_counts

          when: ( 'circrna_discovery' || 'differential_expression' || 'mirna_prediction' in module)

          script:
      		"""
      		python ${projectDir}/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
      		"""
		}

} else{

  single_tool = ciriquant_results.mix(circexplorer2_results, dcc_results, circrna_finder_results, find_circ_results)

  process get_counts_single{

          echo true
          publishDir "$params.outdir/circrna_discovery/count_matrix", mode:'copy'


          input:
            file(bed) from single_tool.collect()
				    val(tool) from params.tool

          output:
            file("circRNA_matrix.txt") into circRNA_counts

          when: ( 'circrna_discovery' || 'differential_expression' || 'mirna_prediction' in module)

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

        when: ('circrna_discovery' || 'mirna_prediction' in module)

        script:
        """
        cp ${projectDir}/bin/unwanted_biotypes.txt ./

        grep -vf unwanted_biotypes.txt $gtf > filt.gtf
        """
}

process get_mature_seq{

        publishDir "$params.outdir/circrna_discovery", mode:'copy', pattern: 'bed12/*.bed'
        publishDir "$params.outdir/circrna_discovery", mode:'copy', pattern: 'fasta/*.fa'

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

        when: ('circrna_discovery' || 'mirna_prediction' in module)

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

        publishDir "$params.outdir/circrna_discovery/annotated", mode: 'copy'

        input:
          file(annotated) from circrna_annotated.collect()

        output:
          file("circrnas_annotated.txt") into annotated_merged

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

	      publishDir "$params.outdir/mirna_prediction/miranda", pattern: "*.miRanda.txt", mode:'copy'

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

    	  publishDir "$params.outdir/mirna_prediction/targetscan", mode:'copy'

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

        publishDir "$params.outdir/mirna_prediction/circos_plots", pattern: "*.pdf", mode: 'copy'
        publishDir "$params.outdir/mirna_prediction/mirna_targets", pattern: "*miRNA_targets.txt", mode: 'copy'

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

        publishDir "$params.outdir/differential_expression/stringtie_quantification", mode:'copy'

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

if(params.phenotype == null){
  exit 1, "[nf-core/circrna] error: parameter '--phenotype' (file for DESeq2) not supplied. Please see '--help' or documentation for help."
}else{
  ch_phenotype = file(params.phenotype)
}

process diff_exp{

        publishDir "$params.outdir/differential_expression", mode:'copy'

	      input:
		      file(gtf_dir) from stringtie_dir.collect()
		      file(circ_matrix) from circrna_matrix_diff_exp
		      file(phenotype) from ch_phenotype

	      output:
		      file("RNA-Seq") into rnaseq_dir
		      file("circRNA") into circrna_dir

        when: 'differential_expression' in module

	      script:
	      """
	      for i in \$(ls -d */); do sample=\${i%"/"}; file=\${sample}.gtf; touch samples.txt; printf "\$sample\t\${i}\${file}\n" >> samples.txt; done

	      prepDE.py -i samples.txt

	      Rscript ${projectDir}/bin/DEA.R gene_count_matrix.csv $phenotype $circ_matrix
	      """
}

(circrna_dir_fetch, circrna_dir_plots) = circrna_dir.into(2)

// obtain the IDs of the differentially expressed circrna

process fetch_de_circ_id{

        input:
          file(circrna_dir) from circrna_dir_fetch

        output:
          file("*.bed") into de_bed_files

        when: 'differential_expression' in module

        script:
        up_reg = "${circrna_dir}/*up_regulated_differential_expression.txt"
      	down_reg = "${circrna_dir}/*down_regulated_differential_expression.txt"
        """
        grep -v "baseMean" $up_reg > up_reg_noheader.txt
        cat $down_reg up_reg_noheader.txt > de_circs.txt

        # make dummy files out of these to place in channel
        awk '{print \$1}' de_circs.txt | grep -v "ID" | while read -r line; do touch \${line}.bed12.bed; done
        """
}

// de circrna dummy bed files in de_bed_files channel, map to get simpleName
ch_de_bed = de_bed_files.flatten().map{ file -> [file.simpleName, file]}

// match incoming bed file channel of all circrnas with DE circ bed ID's
filt_bed_tmp = ch_de_bed.join(bed_diff_exp)

// remove the dummy files (file[1])
filt_bed = filt_bed_tmp.map{file -> [file[0], file[2]]}

// join should keep only the DE keys and apply it to parent, mature files if they are in correct structure.
ch_report = filt_bed.join(parent_diff_exp).join(mature_diff_exp)

// must combine folders here or else process uses once then exits.
ch_DESeq2_dirs = circrna_dir_plots.combine(rnaseq_dir)

process de_plots{

        publishDir "$params.outdir/differential_expression/circrna_expression_plots", pattern:"*.pdf", mode:'copy'

      	input:
      		file(phenotype) from ch_phenotype
      		tuple val(base), file(bed), file(parent_gene), file(mature_length), file(circRNA), file(rnaseq) from ch_report.combine(ch_DESeq2_dirs)

      	output:
      		file("*.pdf") into de_plots
          file("*DESeq2_stats.txt") into de_stats

        when: 'differential_expression' in module

      	script:
      	up_reg = "${circRNA}/*up_regulated_differential_expression.txt"
      	down_reg = "${circRNA}/*down_regulated_differential_expression.txt"
      	circ_counts = "${circRNA}/DESeq2_normalized_counts.txt"
      	gene_counts = "${rnaseq}/DESeq2_normalized_counts.txt"
      	"""
      	# merge upreg, downreg info
        grep -v "baseMean" $up_reg > up_reg_noheader.txt
        cat $down_reg up_reg_noheader.txt > de_circ.txt

      	# Make plots and generate circRNA info
      	Rscript ${projectDir}/bin/circ_report.R de_circ.txt $circ_counts $gene_counts $parent_gene $bed $mature_length $phenotype
      	"""
}

// collect all from previous process
//master_ch = circRNA_plots.collect()
//(test, test1) = circRNA_plots.into(2)
//test.view()
// delete text files in process script, left with only dirs.

process master_report{

        publishDir "$params.outdir/differential_expression/circrna_diff_exp_stats", mode:'copy'

      	input:
      		file(reports) from de_stats.collect()

      	output:
      		file("*circRNAs.txt") into final_out

        when: 'differential_expression' in module

      	script:
      	"""
      	## extract reports
      	for dir in '*/'; do cp \$dir/*_Report.txt .; done

      	# remove header, add manually
      	cat *.txt > merged.txt
      	grep -v "Log2FC" merged.txt > no_headers.txt
      	echo "circRNA_ID Type Mature_Length Parent_Gene Strand Log2FC pvalue Adjusted_pvalue" | tr ' ' '\t' > headers.txt
      	cat headers.txt no_headers.txt > merged_reports.txt

      	Rscript ${projectDir}/bin/annotate_report.R
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
        if ( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable BAM extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

        // output tuple mimicking fromFilePairs if FASTQ provided, else tuple for BAM
        if(bam.matches('NA')){
            [ samples, [read1, read2] ]
        }else{
            [ samples, bam ]
        }

        }

}

// If no input CSV provided, parse input directory containing files.
def retrieve_input_paths(input, glob, type){

      // debug
      println input
      println glob
      println type

      if(type == 'fastq'){

          fastq_files = input + glob
          Channel
              .fromFilePairs(fastq_files)
              .filter{ it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
              .ifEmpty{exit 1, "[nf-core/circrna] error: Your FASTQ files do not have the appropriate extension of either '.fastq.gz', 'fq.gz', fastq' or 'fq'."}
              .map{ row -> [ row[0], [ row[1][0], row[1][1] ]]}
              .ifEmpty{exit 1, "[nf-core/circrna] error: --input was empty - no files supplied"}
              .set{reads_for_csv}

      }else if(type == 'bam'){

          bam_files = input + glob
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
