#!/usr/bin/Rscript

get_args <- function(){

	argp <- arg_parser(
		description="Script to take tool outputs and bring forward circRNAs called by at least 2 tools",
		hide.opts=TRUE)

	argp <- add_argument(
		parser=argp,
		arg="samples",
		short="s",
		help="csv file listing tool output files that are not empty",
		default="samples.csv")
	argv <- parse_args(
		parser=argp,
		argv=commandArgs(trailingOnly=TRUE))


	return(argv)
}

giveError <- function(message){
	cat(paste("\n", message, sep=""))
	quit()
}

usage <- function(){giveError("Usage: consolidate_algorithms.R samples.csv")}

## script body

stage_data <- function(samples){

	inputdata <- list()

	samples <- read.csv(samples, sep="\t", header=F, stringsAsFactors=FALSE)
	names <- gsub(".bed", "", samples$V1)


	inputdata$samples <- samples
	inputdata$names <- names

	return(inputdata)
}


## read in files and make circRNA IDs.

main <- function(inputdata){

	dir.create("tool_id")
	samples <- inputdata$samples
	dflist <- list()

	for(i in 1:nrow(samples)){

		file_handle <- file.path(paste("./", samples$V1[i], sep=""))
		df <- read.table(file_handle, sep="\t", header=F, stringsAsFactors=FALSE)
		dflist[[i]] <- read.table(file_handle, sep="\t", header=F, stringsAsFactors=FALSE)
		circ_id <- with(df, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))
		outfile <- gsub(".bed", ".txt", file_handle)
		write.table(circ_id, paste0("tool_id/", outfile, sep=""), row.names=F, quote=F, sep="\t")
	}

	id_list <- list.files(path="tool_id", pattern=".txt")
	id_list[] <- lapply(id_list, function(x) paste("tool_id/", x, sep=""))

	# create bash command
	vec <- unlist(id_list)
	string <- paste(vec, collapse=" ")

	bash <- paste("cat", string, "| sort | uniq -cd | grep -v \'x\' | awk \'{OFS=\"\t\"; print $2}\' > filteredcirc.txt", sep=" ")
	system(bash)

	filtered <- read.table("filteredcirc.txt", header=F, sep="\t")

	mat <- do.call("rbind", dflist)

	write.table(mat, "combined_counts.bed", sep="\t", row.names = F, col.names = F, quote = F)

}

## error messages, library load
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))

## initiate script
arg <- get_args()

inputdata <- stage_data(arg$samples)

x <- main(inputdata)
