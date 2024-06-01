#!/usr/bin/env Rscript

library(fishpond)

tx_expression <- readRDS('${transcript_rds}')
mi_expression <- read.table('${mirna_expression}', header=TRUE, row.names=1, sep='\\t')
interactions <- read.table('${bindingsites}', sep='\\t')

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version)
    ), 
'versions.yml')
