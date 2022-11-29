//
// Check input files
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet
    ch_phenotype

    main:
    phenotype = params.phenotype && params.module.contains('differential_expression') ? examine_phenotype(ch_phenotype) : Channel.empty()

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    phenotype                                 // Channel: [ pheno ]
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]

        // Conduct check here, cannot figure out how to outside of this scope i.e accessing meta.single end for an if else from channel.
        if(params.tool.split(',').contains('ciriquant')){
            exit 1, "ERROR: Unfortunately CIRIquant does not support single-end reads. Please select one of the other 6 quantification tools."
        }

    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}


def examine_phenotype(pheno){

    Channel
        .fromPath(pheno)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_cols = ['condition']

        if (!row.keySet().containsAll(expected_cols)) exit 1, "[nf-core/circrna] error: 'condition' is not a column name in the phenotype file.\n\nThe primary response variable must be named 'condition', please refer to the usage documentation online"

        def condition  = row.condition.matches('NA') ? 'NA' : row.condition

        if(condition == '') exit 1, "[nf-core/circrna] error: Invalid phenotype file, condition column contains empty cells."
        if(condition.matches('NA')) exit 1, "[nf-core/circrna] error: NA value in phenotype condition column."

        }
        .toList()

        return Channel.value(file(pheno))
}
