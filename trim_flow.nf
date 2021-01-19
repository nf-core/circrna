#!/usr/bin/env nextflow 

// testing trimming logic


params.trimming = null

process trim {

	echo true

	output:
	stdout to out

	when: params.trimming == true

	script:
	"""
	echo "trimming"
	"""
}


