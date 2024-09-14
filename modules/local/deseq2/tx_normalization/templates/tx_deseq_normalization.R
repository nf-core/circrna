#!/usr/bin/env Rscript

library(DESeq2)

raw_counts <- read.table("$counts", sep = "\\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
raw_counts <- raw_counts[ , -2] # drop gene ids
rownames(raw_counts) <- raw_counts\$tx
data <- round(raw_counts[, -1])

samples <- colnames(raw_counts)[-c(1)]


transcript_names <- data.frame(tx = raw_counts\$tx, order = seq_len(nrow(raw_counts)))

# normalize using DeSeq2, Library Size Estimation
meta_data <- data.frame(samples)
row.names(meta_data) <- meta_data\$samples
all(colnames(data) %in% rownames(meta_data))
all(colnames(data) == rownames(meta_data))

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta_data, design = ~ 1)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- DESeq2::counts(dds, normalized = TRUE)

# add tx IDs back to counts table
merged_data <- merge(transcript_names, normalized_counts,
                    by.x = "tx", by.y = "row.names")

merged_data <- merged_data[order(merged_data\$order), ]

norm_data <- subset(merged_data, select = -c(order))

write.table(norm_data, paste0("${meta.id}.normalized_counts.tsv"), quote = FALSE, sep = "\\t", row.names = FALSE)

# TODO: (Can be done later) Add support for Samplesheet so that we can eliminate batch effects


################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-deseq2:', deseq2.version)
    ),
'versions.yml')
