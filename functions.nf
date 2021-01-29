#!/usr/bin/env nextflow

params.input = 'https://raw.githubusercontent.com/BarryDigby/test/dev/test_samples.csv'
params.input_glob = null
params.input_type = null

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


// view it
//ch_input_sample.view()


process foo{

      echo true

      input:
      tuple val(base), file(read) from ch_input_sample

      script:
      """
      echo "$base, $read"
      """
}


// Functions

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
              .view()
              .ifEmpty{exit 1, "[nf-core/circrna] error: --input was empty - no files supplied"}
              .set{reads_for_csv}

      }else if(type == 'bam'){

          bam_files = input + glob
          Channel
              .fromFilePairs(bam_files, size: 1)
              .filter{ it =~/.*.bam/}
              .map{ row -> [row[0], [row[1][0]]]}
              .view()
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
                .view()
}
