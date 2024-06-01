#!/usr/bin/env Rscript

library(fishpond)
suppressMessages(library(SummarizedExperiment))

tx_expression <- readRDS('${transcript_rds}')
mi_expression <- read.table('${mirna_expression}', header=TRUE, row.names=1, sep='\\t')
interactions <- read.table('${bindingsites}', sep='\\t')

result_cols <- c('stat', 'log2FC', 'pvalue', 'locfdr', 'qvalue')

# Iterate rows of interactions
for (i in 1:nrow(interactions)) {
  # Get miRNA and target gene
  miRNA <- interactions[i, 1]
  targets <- unlist(strsplit(interactions[i, 2], ','))

  # TODO: Remove this check after making sure that lowly expressed miRNAs
  # are filtered out before binding site detection
  if (!miRNA %in% rownames(mi_expression)) {
    print(paste('miRNA', miRNA, 'not found'))
    next
  }

  mirna_expression <- mi_expression[miRNA,]
  transcript_expression <- tx_expression[targets,]

  # Add miRNA expression to colData so that it can be used for correlation
  colData(transcript_expression) <- cbind(
    colData(transcript_expression),
    t(mirna_expression[, rownames(colData(transcript_expression))])
  )

  # TODO: Allow setting "spearman"
  result <- rowData(swish(transcript_expression, miRNA, cor = "pearson"))[, result_cols]
  write.table(result, paste0(miRNA, '.tsv'))
}

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
