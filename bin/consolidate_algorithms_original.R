#!/usr/bin/env Rscript

dir = "."
samples = read.table(file.path(dir, 'samples.csv'))
files = file.path(dir, samples[,1])

ciriq_index <- grep("ciriquant", files)
ciriquant <- read.table(files[ciriq_index], header = F, sep ="\t")
ciriquant_id <- with(ciriquant, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))

circrna_finder_index <- grep("circrna_finder", files)
circrna_finder <- read.table(files[circrna_finder_index], header = F, sep="\t")
circrna_finder_id <- with(circrna_finder, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))

circexplorer2_index = grep("circexplorer", files)
circexplorer2 <- read.table(files[circexplorer2_index], header = F, sep ="\t")
circexplorer2_id <- with(circexplorer2, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))

dcc_index = grep("dcc", files)
dcc <- read.table(files[dcc_index], header = F, sep="\t")
dcc_id <- with(dcc, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))

mapsplice_index = grep("mapsplice", files)
mapsplice <- read.table(files[mapsplice_index], header = F, sep="\t")
mapsplice_id <- with(mapsplice, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))

find_circ_index <- grep("find_circ", files)
find_circ <- read.table(files[find_circ_index], header = F, sep="\t")
find_circ_id <- with(find_circ, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))


# write ids to file for bash command
# accept circRNAs in at least 2 tools. 

dir.create("tool_id")
write.table(ciriquant_id, "tool_id/ciriquant.txt", row.names=F, sep="\t", quote=F)
write.table(circexplorer2_id, "tool_id/circexplorer2.txt", row.names=F, sep="\t", quote=F)
write.table(circrna_finder_id, "tool_id/circrna_finder.txt", row.names=F, sep="\t", quote=F)
write.table(dcc_id, "tool_id/dcc.txt", sep="\t", row.names=F, quote=F)
write.table(find_circ_id, "tool_id/find_circ.txt", row.names=F, quote=F, sep="\t")
write.table(mapsplice_id, "tool_id/mapsplice.txt", row.names=F, quote=F, sep="\t")

bash <- 'cat tool_id/circexplorer2.txt tool_id/ciriquant.txt tool_id/circrna_finder.txt tool_id/dcc.txt tool_id/find_circ.txt tool_id/mapsplice.txt | sort | uniq -cd | grep -v \'x\' | awk \'{OFS="\\t"; print $2}\' > filteredcirc.txt'

system(bash)

## read in filtered circRNA IDs
filtered <- read.table("filteredcirc.txt", header=F, sep="\t")

mat <- rbind(ciriquant, circexplorer2, circrna_finder, find_circ, dcc, mapsplice)
mat$id <- with(mat, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))
mat <- mat[which(mat$id %in% filtered$V1),]
mat1 <- mat[order(mat[,6], -abs(mat[,5])),]
mat1 <- mat1[!duplicated(mat1$id),]

mat1 <- mat1[,1:5]

write.table(mat1, "combined_counts.bed", sep="\t", row.names = F, col.names = F, quote = F)
