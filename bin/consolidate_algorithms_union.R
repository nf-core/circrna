#!/usr/bin/Rscript

get_args <- function(){

    argp <- arg_parser(
        description="Script to take union of quant tools. If duplicate IDs, takes highest count.",
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

##
## Strategy for union:
## 1. read in non empty bed files using samples.csv loop
## 2. append each tool output to dflist
## 3. convert to master DF (rbind)
## 4. make circ ID column
## 5. sort by ID and abs count value
## 6. remove duplicate IDs (keep highest count due to sort)
## 7. remove ID column, write to file.

stage_data <- function(samples){

    inputdata <- list()
    samples <- read.csv(samples, sep="\t", header=F, stringsAsFactors=FALSE)
    inputdata$samples <- samples
    return(inputdata)
}
main <- function(inputdata){

    samples <- inputdata$samples
    dflist <- list()

    for(i in 1:nrow(samples)){

        file_handle <- file.path(paste("./", samples$V1[i], sep=""))
        df <- read.table(file_handle, sep="\t", header=F, stringsAsFactors=FALSE)
        dflist[[i]] <- read.table(file_handle, sep="\t", header=F, stringsAsFactors=FALSE)

    }

    master <- do.call("rbind", dflist)
    master$ID <- with(master, paste0(V1, sep="-", V2, sep=":", V3, sep="-", V4))

    master_sort <- master[order(master$ID, -abs(master$V5) ), ]
    master_sort <- master_sort[ !duplicated(master_sort$ID), ]

    # drop ID column and write to file
    master_union <- subset(master_sort, select=-c(ID))

    write.table(master_union, "combined_counts.bed", sep="\t", row.names = F, col.names = F, quote = F)

}

## error messages, library load
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))

## initiate script
arg <- get_args()

inputdata <- stage_data(arg$samples)

x <- main(inputdata)
