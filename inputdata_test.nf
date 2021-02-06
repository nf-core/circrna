#!/usr/bin/env nextflow

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
    inputSample = retrieve_input_paths(params.input, params.input_type)
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

        label 'py3'

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
                out1=${base}_R1.trim.fq.gz \
                out2=${base}_R2.trim.fq.gz \
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

                label 'py3'

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

              	label 'py3'

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


// Stage Aligner read channels
(circexplorer2_reads, find_circ_reads, ciriquant_reads, mapsplice_reads, uroborus_reads, circrna_finder_reads, dcc_reads, dcc_reads_mate1, dcc_reads_mate2, hisat2_reads) = aligner_reads.into(10)
dcc_reads.view()


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
