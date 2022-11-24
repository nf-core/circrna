#!/usr/bin/Rscript

get_args <- function(){

    argp <- arg_parser(
        description="Script to take tool outputs and bring forward circRNAs called by at least n tools",
        hide.opts=TRUE)

    argp <- add_argument(
        parser=argp,
        arg="samples",
        short="s",
        help="csv file listing tool output files that are not empty",
        default="samples.csv")

    argp <- add_argument(
        parser=argp,
        arg="n_tools",
        short="n",
        help="number of tools circrna must be called by")

    argp <- add_argument(
        parser=argp,
        arg="duplicates_fun",
        short="d",
        help="how to handle counts for duplicate circRNAs [mean, max]")

    argv <- parse_args(
        parser=argp,
        argv=commandArgs(trailingOnly=TRUE))

    return(argv)
}

giveError <- function(message){
    cat(paste("\n", message, sep=""))
    quit()
}

usage <- function(){giveError("Usage: consolidate_algorithms.R samples.csv ${params.n_tools}")}

## script body

stage_data <- function(samples, n_tools, duplicates_fun){

    inputdata <- list()

    samples <- read.csv(samples, sep="\t", header=F, stringsAsFactors=FALSE)
    inputdata$samples <- samples
    inputdata$n_tools <- n_tools
    inputdata$duplicates_fun <- duplicates_fun

        return(inputdata)
}


## read in files and make circRNA IDs.

main <- function(inputdata){

    samples <- inputdata$samples
    n_tools <- inputdata$n_tools
    duplicates_fun <- inputdata$duplicates_fun

    ## intit lists
    dflist <- list()
    idlist <- list()

    # loop over samples, append to counts mat and IDs to respective lists
    for(i in 1:nrow(samples)){
        file_handle <- file.path(paste("./", samples$V1[i], sep=""))
        df <- read.table(file_handle, sep="\t", header=F, stringsAsFactors=FALSE)
        dflist[[i]] <- read.table(file_handle, sep="\t", header=F, stringsAsFactors=FALSE)
        idlist[[i]] <- with(df, paste0(V1, sep="-", V2, sep=":", V3, sep="-", V4))
    }

    # place all ids in a vector
    vec <- unlist(idlist)

    # make table to get ID ocurrence count
    tab <- table(vec)

    # now filter by ocurrence i.e must have been in "n" tools
    filt_id <- names(tab)[tab >= n_tools]

    # make matser df, append ID
    mat <- do.call("rbind", dflist)
    mat$ID <- with(mat, paste0(V1, sep="-", V2, sep=":", V3, sep="-", V4))

    # extract circrnas called by n tools
    mat <- subset(mat, mat$ID %in% filt_id)
    print(head(mat))
    # handle duplicate circRNA counts:
    if(duplicates_fun == "max"){

        # Take the max value for duplicate circRNAs
        mat <- mat[order(mat$ID, -abs(mat$V5)),]
        mat <- mat[!duplicated(mat$ID),]
        mat <- subset(mat, select=-c(ID))
    }else{

        # Remove duplicates like above to generate row indices
        tmp <- mat[order(mat$ID, -abs(mat$V5)),]
        tmp <- tmp[!duplicated(tmp$ID),]
        # Now Take average values, append back to tmp.
        mat <- setNames(aggregate(mat$V5, list(mat$ID), FUN=mean), c("ID", "V5"))
        stopifnot(nrow(tmp)==nrow(mat))
        tmp <- subset(tmp, select=c(V1,V2,V3,V4))
        mat <- cbind(tmp,mat)
        # round up to nearest integer
        mat$V5 <- ceiling(mat$V5)
        mat <- subset(mat, select=-c(ID))
    }

    write.table(mat, "combined_counts.bed", sep="\t", row.names = F, col.names = F, quote = F)

}

## error messages, library load
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))

## initiate script
arg <- get_args()

inputdata <- stage_data(arg$samples, arg$n_tools, arg$duplicates_fun)

x <- main(inputdata)
