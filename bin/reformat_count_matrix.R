#!/usr/bin/env Rscript

## Author: Barry Digby
## License: MIT

library(dplyr)
mat <- read.table("circRNA_matrix.txt", sep="\t", header=T, stringsAsFactors=F, check.names=F)
mat$ID <- with(mat, paste0(Chr, sep=":", Start, sep="-", Stop, sep=":", Strand))
mat <- mat[,-c(1:4)]
mat1 <- mat %>% select(ID, everything())
ID <- as.data.frame(mat1$ID)
mat <- as.data.frame(subset(mat1, select=-c(ID)))
mat <- mat[, order(names(mat))]
mat1 <- cbind(ID, mat)
colnames(mat1)[1] <- "ID"
write.table(mat1, "count_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F)
