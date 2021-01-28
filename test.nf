#!/usr/bin/env nextflow

params.index = 'test_samples.csv'
params.flow = null

Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleID, file(row.read1), file(row.read2)) }
    .into{ samples_ch1; samples_ch }

samples_ch1.view()


if(params.flow == 'A'){

	println "Echo"

}else if(params.flow == 'B'){

	println "Exho"

}else if(params.flow == 'C'){

	println "try"
}


