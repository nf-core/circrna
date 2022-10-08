#!/usr/bin/Rscript

circ_mat = read.table("circRNA_matrix.txt", header=T, sep="\t", check.names = FALSE, stringsAsFactors = F, row.names = "ID")
gene_mat = read.table("gene_count_matrix.csv", sep=",", header=T, row.names="gene_id", stringsAsFactors = F)
map = read.table("circrna_host-gene.txt", header = F, sep="\t", stringsAsFactors = F)

# some circrnas do not have a host gene.
map <- na.omit(map)
colnames(map) <- c("circrna", "gene")

# resolve multiple host genes by creating a new 'map' dataframe
# to enforce 1-1 mapping.

new_circ = c()
new_gene = c()

for(i in 1:nrow(map)){
  row <- map[i,]
  circ <- row$circrna
  gene <- row$gene

  # multiple host genes?
  multiple_genes <- unlist(strsplit(gene, ","))
  #print(multiple_genes)
  if(length(multiple_genes) > 1){
    for(gene in multiple_genes){
      new_circ <- c(new_circ, circ)
      new_gene <- c(new_gene, gene)
    }
  }
  new_circ <- c(new_circ, circ)
  new_gene <- c(new_gene, gene)
}

new_map <- data.frame(new_circ, new_gene)

#create circTest dataframe. really odd formatting..!
new_circ_mat <- circ_mat[c(new_map$new_circ),]

# unique names for duplicated circrna ids, resolve below
Chr <- c()
Start <- c()
End <- c()
for(i in 1:nrow(new_circ_mat)){
  row <- new_circ_mat[i,]
  print(rownames(row))
  chr_ <- unlist(strsplit(rownames(row), ':'))[1]
  coords <- unlist(strsplit(rownames(row), ':'))[2]
  start_ <- unlist(strsplit(coords, '-'))[1]
  end_ <- unlist(strsplit(coords, '-'))[2]

  Chr <- c(Chr, chr_)
  Start <- c(Start, start_)
  End <- c(End, end_)

}

Gene <- new_map$new_gene
circ_csv <- data.frame(Chr, Start, End, Gene, new_circ_mat)
# DROP the rownames, do not write to file.

gene_csv <- gene_mat[c(new_map$new_gene),]
gene_csv <- data.frame(Chr, Start, End, Gene, gene_csv)

write.csv(circ_csv, "circ.csv", quote=F, row.names = FALSE)
write.csv(gene_csv, "linear.csv", quote=F, row.names = FALSE)
