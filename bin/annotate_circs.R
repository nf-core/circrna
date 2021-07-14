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

usage <- function(){giveError("USAGE: annotate_circs.R parent_gene bed mature_length")}

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

	write.table(out_df, paste(file_name, "annotated.txt", sep="_"), quote=F, sep="\t", row.names=F)
}

# Packages + Error traceback (really handy for explicit error tracing)
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("circlize"))

arg <- get_args()
inputdata <- stage_data(arg$parent_gene, arg$bed, arg$mature_len)
z <- singular_report(inputdata)
