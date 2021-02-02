#!/usr/bin/Rscript

get_args <- function(){

	argp <- arg_parser(
		description="output annotated circrna incorporating mature len, parent gene info, circrna type info",
		hide.opts=TRUE)

	argp <- add_argument(
		parser=argp,
		arg="parent_gene",
		short="pg",
		help="circRNA parent gene file")

	argp <- add_argument(
		parser=argp,
		arg="bed",
		short="b",
		help="bed12 file of circRNA")

	argp <- add_argument(
		parser=argp,
		arg="mature_len",
		short="l",
		help="length of mature circRNA")

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

stage_data <- function(parent_gene, bed, mature_len){

	inputdata <- list()

	parent_tmp <- read.table(parent_gene, sep="\t", row.names=1, header=F, stringsAsFactors=F)
	colnames(parent_tmp) <- "gene"
	parent <- parent_tmp$gene
	bed <- read.table(bed, sep="\t", header=F, stringsAsFactors=F)
	colnames(bed) <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "type", "ExonCount", "ExonSizes", "ExonStart")
	mature_tmp <- read.table(mature_len)
	mature <- mature_tmp$V1

	inputdata$parent <- parent
	inputdata$bed <- bed
	inputdata$mature <- mature

	return(inputdata)
}

data_summary <- function(data, varname, groupnames){
  	require(plyr)
  	summary_func <- function(x, col){
    		c(mean = mean(x[[col]], na.rm=TRUE), sd = sd(x[[col]], na.rm=TRUE))
		}
  	data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
  	data_sum <- rename(data_sum, c("mean" = varname))

	return(data_sum)
}


singular_report <- function(inputdata){

	chr <- inputdata$bed$chr
	start <- inputdata$bed$start
	end <- inputdata$bed$end
	coords <- paste(start, end, sep="-")
	circ_id <- paste(chr, coords, sep=":")
	strand <- inputdata$bed$strand
	gene <- inputdata$parent
	mature_len <- inputdata$mature
	file_name <- inputdata$bed$name
	type <- inputdata$bed$type

	vec <- c(circ_id, type, mature_len, gene, strand)
	mat <- matrix(vec, ncol=5)
	out_df <- as.data.frame(mat, stringsAsFactors=FALSE)

	colnames(out_df) <- c("circRNA_ID", "Type", "Mature_Length", "Parent_Gene", "Strand")

	write.table(out_df, file.path(file_name, paste(file_name, "annotated.txt", sep="_")), quote=F, sep="\t", row.names=F)
}

# Packages + Error traceback (really handy for explicit error tracing)
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("circlize"))

arg <- get_args()
inputdata <- stage_data(arg$de_circ, arg$circ_counts, arg$gene_counts, arg$parent_gene, arg$bed, arg$miranda, arg$targetscan, arg$mature_len, arg$phenotype, arg$circlize_exons)
dir.create(inputdata$bed$name)
x <- prep_plots(inputdata)
y <- miRNAs(inputdata)
z <- singular_report(inputdata)

#head(inputdata$de)
#head(inputdata$circ)
#head(inputdata$gene)
#head(inputdata$parent)
#head(inputdata$bed)
#head(inputdata$miranda)
#head(inputdata$targetscan)
#head(inputdata$mature)
#head(inputdata$pheno)
