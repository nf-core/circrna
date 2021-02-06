#!/usr/bin/Rscript

library(biomaRt)

x <- read.table("merged_reports.txt", sep="\t", header=T)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes=c("entrezgene_description", "hgnc_symbol"), mart=mart, useCache=FALSE)
colnames(results) <- c("Description", "Parent_Gene")

# allow for missing descriptions! Crucial to capture lncRNAs
x <- merge(x, results, by="Parent_Gene", all.x=T)
x <- subset(x, select=c(circRNA_ID, Type, Mature_Length, Parent_Gene, Strand, Log2FC, pvalue, Adjusted_pvalue, Description))

up <- subset(x, x$Log2FC > 0)
down <- subset(x, x$Log2FC < 0)

up <- up[order(abs(up$Log2FC), decreasing=T),]
down <- down[order(abs(down$Log2FC), decreasing=T),]

write.table(up, "Up_Regulated_circRNAs.txt", sep="\t", quote=F, row.names=F)
write.table(down, "Down_Regulated_circRNAs.txt", sep="\t", quote=F, row.names=F)
write.table(x, "DE_circRNAs.txt", sep="\t", quote=F, row.names=F)
