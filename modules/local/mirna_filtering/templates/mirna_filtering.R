#!/usr/bin/env Rscript

expression_norm <- read.table("$normalized_counts",
    sep = "\\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

samples <- colnames(expression_norm)[-c(1)]

# filter data: counts > 5 in at least 20% of samples
if (length(samples) < 5) {
    stop("Cannot perform filtering on less than 5 samples")
}

sample_nr_cutoff <- ceiling($mirna_min_sample_percentage * length(samples))
rows_to_keep <- c()

for (i in seq_len(nrow(expression_norm))) {
  mirna_per_sample <- 0
  for (j in 5:ncol(expression_norm)) {
    if (expression_norm[i, j] >= $mirna_min_reads) {
      mirna_per_sample <- mirna_per_sample + 1
    }
  }
  if (mirna_per_sample >= sample_nr_cutoff) {
    rows_to_keep <- append(rows_to_keep, i)
  }
}

filtered_data <- expression_norm[rows_to_keep, ]

write.table(filtered_data, paste0("${meta.id}.normalized_counts_filtered.tsv"),
            quote = FALSE, sep = "\\t",
            row.names = FALSE)

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
