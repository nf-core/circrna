#!/usr/bin/env Rscript

circ <- read.table("${circ_counts}", header=T, sep="\\t", check.names = FALSE)
gene <- read.table("${gene_counts}", header=T, sep="\\t", check.names = FALSE, row.names = 1)

gene <- gene[circ\$gene_id, ]

rownames(circ) <- circ\$tx
rownames(gene) <- rownames(circ)
circ\$tx <- NULL

if (nrow(circ) != nrow(gene)) {
    stop("Number of rows in circ and gene counts do not match")
}

if (nrow(circ) > 0) {
    write.table(circ, "${prefix}_circs.tsv", sep="\\t", quote=F, row.names=T)
    write.table(gene, "${prefix}_genes.tsv", sep="\\t", quote=F, row.names=T)
}

writeLines(c(
    '"${task.process}":',
    paste0('    r-base: "', getRversion(), '"')
), 'versions.yml')
