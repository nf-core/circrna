#!/usr/bin/env Rscript --vanilla

# Script for importing and processing transcript-level quantifications.
# Written by Lorena Pantano, later modified by Jonathan Manning, and released
# under the MIT license.

# Loading required libraries
library(tximeta)

################################################
################################################
## Main script starts here                    ##
################################################
################################################

# Define pattern for file names based on quantification type
pattern <- ifelse('$quant_type' == "kallisto",
                ifelse(length(list.files('quants', pattern = "abundance.h5", recursive = T, full.names = T)) != 0,
                    "abundance.h5",
                    "abundance.tsv"),
                "quant.sf")

fns <- list.files('quants', pattern = pattern, recursive = T, full.names = T)
names <- basename(dirname(fns))
names(fns) <- names

coldata <- data.frame(files = fns, names = names)
rownames(coldata) <- coldata[["names"]]

# Import transcript-level quantifications
se <- tximeta(coldata, type = '$quant_type', txOut = TRUE)

# Save summarized experiment to file
saveRDS(se, file = paste0('$prefix', '.rds'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste("R_sessionInfo.log", sep = '.'))
citation("tximeta")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
tximeta.version <- as.character(packageVersion('tximeta'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-tximeta:', tximeta.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
