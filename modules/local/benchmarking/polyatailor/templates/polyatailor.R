#!/usr/bin/env Rscript
library(PolyAtailor)
library(Biostrings)
library(yaml) # To write the YAML file
library(Rsamtools) # Ensure Rsamtools is loaded for scanBam
library(GenomicAlignments)
library(dplyr)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("strsplit", "base")

fastq_string <- '$fastq'
fastq_files <- unlist(strsplit(fastq_string, " "))
sample_names <- gsub(".fq.gz\$", "", fastq_files)

result_1 <- tailScan(
  fastq = fastq_files[1],
  mcans = 3,                      
  findUmi = FALSE,                
  lumi = 0,                       
  adapterSeq = "",                
  anchorSeq = "",                 
  resultpath = "./",              
  samplename = sample_names[1],
  tailAnchorLen = 8,              
  minTailLen = 5,                 
  realTailLen = 15,               
  maxNtail = 2,                   
  mapping = FALSE,                
  mapinfo = NULL,                 
  findTailType = 'A'              
)

write.table(result_1, file = paste0(sample_names[1], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)

result_2 <- tailScan(
  fastq = fastq_files[2],
  mcans = 3,                      
  findUmi = FALSE,                
  lumi = 0,                       
  adapterSeq = "",                
  anchorSeq = "",                 
  resultpath = "./",              
  samplename = sample_names[2],
  tailAnchorLen = 8,              
  minTailLen = 5,                 
  realTailLen = 15,               
  maxNtail = 2,                   
  mapping = FALSE,                
  mapinfo = NULL,                 
  findTailType = 'A'              
)
write.table(result_2, file = paste0(sample_names[2], ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)

result_tailMap <- tailMap(
  bamfile = '$bam',
  mcans = 3,
  minTailLen = 5,
  findUmi = FALSE,
  maxNtail = 2,
  mapping = TRUE,
  longRead = FALSE
)
write.csv(result_tailMap, "map.csv", row.names = FALSE)



version_info <- list(
  PolyAtailor_version = as.character(packageVersion("PolyAtailor")),
  Biostrings_version = as.character(packageVersion("Biostrings")),
  run_time = Sys.time()
)
write_yaml(version_info, "versions.yml")
