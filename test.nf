#!/usr/bin/env nextflow

params.index = 'test_samples.csv'

Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleID, row.read1, row.read2) }
    .into{ samples_ch1; samples_ch }

samples_ch1.view()
