#!/usr/bin/Rscript

get_args <- function(){

  argp <- arg_parser(
    description="output circrna-mirna circos plot, combined miranda + targetscan outputs to file",
    hide.opts=TRUE)

  argp <- add_argument(
    parser=argp,
    arg="bed",
    short="b",
    help="bed12 file of circRNA")

  argp <- add_argument(
    parser=argp,
    arg="miranda",
    short="m",
    help="miRanda output for circRNA")

  argp <- add_argument(
    parser=argp,
    arg="targetscan",
    short="t",
    help="targetscan output for circRNA")

  argp <- add_argument(
    parser=argp,
    arg="circlize_exons",
    short="z",
    help="output from bash script preparing circlize exons")

  argp <- add_argument(
    parser=argp,
    arg="species_id",
    short="s",
    help="grab this from miRanda file using bash and supply to script")

  argv <- parse_args(
    parser=argp,
    argv = commandArgs(trailingOnly = TRUE))

  return(argv)
}

giveError <- function(message){
  cat(paste("\n", message, sep=""))
  quit()
}

usage <- function(){giveError("USAGE: circ_report.R bed miranda targetscan mature_len circlize_exons.txt")}

stage_data <- function(bed, miranda, targetscan, circlize_exons, species){

  inputdata <- list()
  bed <- read.table(bed, sep="\t", header=F, stringsAsFactors=F)
  colnames(bed) <- c("chr", "start", "end", "name", "bsj_reads", "strand", "thickStart", "thickEnd", "ExonCount", "ExonSizes", "ExonStart", "type", "genes", "transcripts", "mature_length")
  miranda <- read.table(miranda, sep="\t", header=T, stringsAsFactors=F)
  targetscan <- read.table(targetscan, sep="\t", header=T, stringsAsFactors=F)
  circlize_exons <- read.table(circlize_exons, sep="\t")
  species <- as.character(species)

  inputdata$bed <- bed
  inputdata$miranda <- miranda
  inputdata$targetscan <- targetscan
  inputdata$circlize_exons <- circlize_exons
  inputdata$species <- species
  return(inputdata)
}

miRNAs <- function(inputdata){

  outdir <- inputdata$bed$name

  miranda <- inputdata$miranda
  targetscan <- inputdata$targetscan
  x <- inputdata$circlize_exons
  circlize_exons <- inputdata$circlize_exons

  miranda$miRNA <- gsub(inputdata$species, "", miranda$miRNA)
  colnames(targetscan)[2] <- "miRNA"

  mir_df <- merge(miranda, targetscan, by="miRNA")

  ## Filtering Step:
  ## Removes miRNA <= specified MFE filter value
  ## Removes duplicate miRNAs sharing the same bind site (duplicate miRNA IDs)
  ## Keeps the miRNA with the higher "Score"

  ## Add a check here, if abs(MFE) is greater than largest energy value, discard filtering method
  #if(abs(as.numeric(mfe)) > min(mir_df$Energy_KcalMol)){
    mir_df <- mir_df[order(mir_df$MSA_start, -abs(mir_df$Score)),]
    mir_df <- mir_df[!duplicated(mir_df$MSA_start),]
    miRs <- subset(mir_df, select=c(miRNA, MSA_start, MSA_end))
    colnames(miRs) <- c("miRNA", "Start", "End")
  #}else{
   # mir_df <- subset(mir_df, mir_df$Energy_KcalMol <= as.numeric(mfe))
  #  mir_df <- mir_df[order(mir_df$MSA_start, -abs(mir_df$Score)),]
   # mir_df <- mir_df[!duplicated(mir_df$MSA_start),]
  #  miRs <- subset(mir_df, select=c(miRNA, MSA_start, MSA_end))
   # colnames(miRs) <- c("miRNA", "Start", "End")
  #}

  ## Calculate where the miRNAs fall in the context of the exons of the circRNA
  ## Add ifelse to allow for "empty" exons with no miRs.
  if(nrow(x)==1){
    exon_1 <- miRs[which(miRs$Start < x$V3[1]),]
    exon_1$value = 1
    exon_1$chr <- "exon1"
    exon_1 <- exon_1[,c(5,2,3,4,1)]
    circlize_mirs <- exon_1
  }else{
    circlize_mirs <- data.frame(ncol(5))
    for(n in 1:nrow(x)){
      if(n==1){
        exon_1 <- miRs[which(miRs$Start < x$V3[n]),]
        if(dim(exon_1)[1]==0){
          rm(exon_1)
        }else{
          exon_1$value = 1
          exon_1$chr <- "exon1"
          exon_1 <- exon_1[,c(5,2,3,4,1)]
          circlize_mirs <- exon_1
        }
      }else{
        exon <- miRs[which(miRs$Start >=  sum(x$V3[0:(n-1)]) & miRs$Start < sum(x$V3[0:n])),]
        if(dim(exon)[1]==0){
          rm(exon)
        }else{
          exon$value = 1
          exon$chr <- paste("exon", n, sep="")
          subtract_me <- sum(x$V3[0:(n-1)])
          exon$Start <- (exon$Start) - subtract_me
          exon$End <- (exon$End) - subtract_me
          exon <- exon[,c(5,2,3,4,1)]
          circlize_mirs <- rbind(circlize_mirs, exon)
        }
      }
    }
  }

  write_mirs <- circlize_mirs
  write_mirs <- subset(write_mirs, select=-c(value, Start, End))
  write_mirs <- cbind(write_mirs, mir_df)
  write_mirs <- subset(write_mirs, select=c(miRNA, Score, Energy_KcalMol, MSA_start, MSA_end, Site_type))
  colnames(write_mirs) <- c("miRNA", "Score", "Energy_KcalMol", "Start", "End", "Site_type")

  write.table(write_mirs, paste(outdir, "miRNA_targets.txt", sep="_"), sep="\t", row.names=F, quote=F)
  make_circos_plot(inputdata, circlize_exons, circlize_mirs)
}

make_circos_plot <- function(inputdata, circlize_exons, circlize_mirs){

  circ_id <- inputdata$bed$name
  col_text <- "grey25"

  pdf(paste(circ_id, "miRNA_Plot.pdf", sep="_"))
  circos.initialize(factors=circlize_exons$V1, xlim=matrix(c(circlize_exons$V2,circlize_exons$V3), ncol=2))
  circos.genomicLabels(circlize_mirs,labels.column = 5, side = "outside", niceFacing = TRUE, cex = 0.8)
  circos.track(ylim=c(0,0.5), panel.fun=function(x,y){
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),chr)
  })
  brk <- seq.int(0,100000,25)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(h="top",major.at=brk,labels=round(brk/1,1),labels.cex=0.4,
                col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")},bg.border=F)
  dev.off()
}


# Packages + Error traceback (really handy for explicit error tracing)
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("circlize"))

arg <- get_args()
inputdata <- stage_data(arg$bed, arg$miranda, arg$targetscan, arg$circlize_exons, arg$species)
y <- miRNAs(inputdata)

#head(inputdata$de)
#head(inputdata$circ)
#head(inputdata$gene)
#head(inputdata$parent)
#head(inputdata$bed)
#head(inputdata$miranda)
#head(inputdata$targetscan)
#head(inputdata$mature)
#head(inputdata$pheno)
