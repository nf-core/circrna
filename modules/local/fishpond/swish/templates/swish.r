#!/usr/bin/env Rscript --vanilla

library(fishpond)
library(SummarizedExperiment)

paths <- c('${experiments.join("\', \'")}')
experiments <- lapply(paths, readRDS)

phenotype <- read.csv('${phenotype}', stringsAsFactors = FALSE)

se_assays <- list()

for (se in experiments) {
  assays <- assays(se)
  # Iterate over named list of assays
  for (assay_name in names(assays)) {
    assay <- assays[[assay_name]]
    
    # Add assay to se_assays for its name
    if (is.null(se_assays[[assay_name]])) {
      se_assays[[assay_name]] <- assay
    } else {
      se_assays[[assay_name]] <- cbind(se_assays[[assay_name]], assay)
    }
  }
}

se_cbind <- do.call(SummarizedExperiment::cbind, experiments)
se <- SummarizedExperiment(assays = se_assays, colData = colData(se_cbind), rowData = rowData(se_cbind))

# Join phenotype data
colData(se) <- merge(colData(se), phenotype, by.x="names", by.y="Sample_ID")
colData(se)\$condition <- as.factor(colData(se)\$condition)

se <- scaleInfReps(se)
se <- labelKeep(se)
se <- swish(se, x="condition")

saveRDS(se, '${meta.id}.rds')

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-fishpond:', packageVersion('fishpond'))
    ),
'versions.yml')