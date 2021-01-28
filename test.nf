#!/usr/bin/env nextflow

params.index = 'test_samples.csv'
params.flow = null

Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleID, [file(row.read1), file(row.read2)]) }
    .set{ samples_ch }

(samples_view, fastqc_reads, raw_reads) = samples_ch.into(3)
samples_view.view()

if(params.flow == 'A'){

	println "Echo"

}else if(params.flow == 'B'){

	println "Exho"

}else if(params.flow == 'C'){

	println "try"
}

process fastqc{

	publishDir ".", mode:'copy'

	input:
	tuple val(base), file(reads) from fastqc_reads

	output:
	file("*{.html,.zip}") into fastqc_out

	script:
	"""
	fastqc -q $reads
	"""
}

aligner_reads = raw_reads

aligner_reads.view()
