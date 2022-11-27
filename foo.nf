#!/usr/bin/env nextflow 

a = Channel.empty()

b = Channel.fromPath('main.nf')

a.mix(b).view()
