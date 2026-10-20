#!/usr/bin/env Rscript

library(fishpond)
suppressMessages(library(SummarizedExperiment))

tx_expression <- readRDS('${transcript_rds}')
mi_expression <- read.table('${mirna_expression}', header=TRUE, row.names=1, sep='\\t', check.names=FALSE)
interactions <- read.table('${bindingsites}', sep='\\t')

missing_samples <- setdiff(rownames(colData(tx_expression)), colnames(mi_expression))
if (length(missing_samples) > 0) {
    stop(paste0(
        'RNA-seq samples missing from the miRNA expression table: ', paste(missing_samples, collapse = ', '),
        '. Available miRNA samples: ', paste(colnames(mi_expression), collapse = ', ')
    ))
}

tx_expression <- scaleInfReps(tx_expression)
tx_expression <- labelKeep(tx_expression) # Here one can perform custom filtering

if (!any(mcols(tx_expression)\$keep)) {
    stop('No transcripts left after filtering')
}

result_cols <- c('stat', 'log2FC', 'pvalue', 'locfdr', 'qvalue')

# Iterate rows of interactions
for (i in 1:nrow(interactions)) {
    # Get miRNA and target gene
    miRNA <- interactions[i, 1]
    targets <- unlist(strsplit(interactions[i, 2], ','))

    mirna_expression <- mi_expression[miRNA,]
    transcript_expression <- tx_expression[targets,]

    if (!any(mcols(transcript_expression)\$keep)) {
        print(paste('No transcripts left after filtering for miRNA', miRNA))
        next
    }

    # Add miRNA expression to colData so that it can be used for correlation
    colData(transcript_expression) <- cbind(
        colData(transcript_expression),
        t(mirna_expression[, rownames(colData(transcript_expression))])
    )

    result <- rowData(swish(transcript_expression, miRNA, cor = "${params.mirna_correlation}"))[, result_cols]
    result <- result[complete.cases(result), ]
    write.table(result, paste0(miRNA, '.tsv'), sep = '\\t')
}

writeLines(c(
    '"${task.process}":',
    paste0('    r-base: "', getRversion(), '"'),
    paste0('    bioconductor-fishpond: "', packageVersion("fishpond"), '"')
), 'versions.yml')
