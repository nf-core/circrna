#!/usr/bin/Rscript

get_args <- function(){

	argp <- arg_parser(
		description="this script is designed to produce plots of DE circRNAs detected in samples as part of the circRNA nextflow script",
		hide.opts=TRUE)

	argp <- add_argument(
		parser=argp,
		arg="de_circ",
		short="d",
		help="DESeq2 Results for DE circs")

	argp <- add_argument(
    parser=argp,
    arg="circ_counts",
		short="c",
    help="Normalized circRNA counts")

 	argp <- add_argument(
    parser=argp,
    arg="gene_counts",
		short="g",
    help="Normalized gene counts")

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

	argp <- add_argument(
		parser=argp,
		arg="phenotype",
		short="p",
		help="Phenotype / colData / samples file used previously for DESeq2")

    	argv <- parse_args(
            	parser=argp,
            	argv = commandArgs(trailingOnly = TRUE))

    return(argv)
    }

giveError <- function(message){
    cat(paste("\n", message, sep=""))
    quit()
    }

usage <- function(){giveError("USAGE: circ_report.R de_circ circrna_counts gene_counts bed mature_len phenotype")}

stage_data <- function(de_circ, circ_counts, gene_counts, parent_gene, bed, mature_len, phenotype){

	inputdata <- list()

	de <- read.table(de_circ, sep="\t", row.names="ID", header=T, stringsAsFactors=F)
	circ <- read.table(circ_counts, sep="\t", row.names="ID", header=T)
	gene <- read.table(gene_counts, sep="\t", row.names="ID", header=T, check.names=F)
	parent_tmp <- read.table(parent_gene, sep="\t", row.names=1, header=F, stringsAsFactors=F)
	colnames(parent_tmp) <- "gene"
	parent <- parent_tmp$gene
	bed <- read.table(bed, sep="\t", header=F, stringsAsFactors=F)
	colnames(bed) <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "type", "ExonCount", "ExonSizes", "ExonStart")
	mature_tmp <- read.table(mature_len)
	mature <- mature_tmp$V1
	pheno <- read.table(phenotype, sep="\t", header=T, row.names=1)

	inputdata$de <- de
	inputdata$circ <- circ
	inputdata$gene <- gene
	inputdata$parent <- parent
	inputdata$bed <- bed
	inputdata$mature <- mature
	inputdata$pheno <- pheno

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

	de_df <- inputdata$de
	circ_de_info <- de_df[which(rownames(de_df) %in% file_name),]
	numcols <- c("log2FoldChange", "pvalue", "padj")
  circ_de_info[numcols] <- sapply(circ_de_info[numcols],as.numeric)
	log2fc <- round(circ_de_info$log2FoldChange, digits=3)
	pval <- signif(circ_de_info$pvalue, digits=3)
	adj.pval <- signif(circ_de_info$padj, digits=3)

	vec <- c(circ_id, type, mature_len, gene, strand, log2fc, pval, adj.pval)
	mat <- matrix(vec, ncol=8)
	out_df <- as.data.frame(mat, stringsAsFactors=FALSE)

	colnames(out_df) <- c("circRNA_ID", "Type", "Mature_Length", "Parent_Gene", "Strand", "Log2FC", "pvalue", "Adjusted_pvalue")

	write.table(out_df, paste(file_name, "DESeq2_stats.txt", sep="_"), quote=F, sep="\t", row.names=F)
}

prep_plots <- function(inputdata){

	circ_key <- inputdata$bed$name
	rna_key <- inputdata$parent

	circ <- as.data.frame(inputdata$circ)
	rna <- as.data.frame(inputdata$gene)
	pheno <- as.data.frame(inputdata$pheno)

	circ_df <- circ[which(rownames(circ) %in% circ_key),]
	rna_df <- rna[which(rownames(rna) %in% rna_key),]

	circ_df <- t(circ_df)
	rna_df <- t(rna_df)

	# it is possible that the parent gene is not in the RNA matrix
	# the df is transposed, check if columns are ifEmpty

  # boxplots dont need rna and satisfy process output (*.pdf)
	circ_df <- cbind(circ_df, pheno)
	circ_df$type <- "circRNA"

	circ_names <- colnames(circ_df)
	circ_names[1] <- "count"
	colnames(circ_df) <- circ_names

	make_boxplot(circ_df, inputdata)

	if(dim(rna_df)[2] == 0){

			print("Parent Gene not present in RNA-Seq matrix, skipping plots")

		}else{

		rna_df <- cbind(rna_df, pheno)

		rna_df$type <- "linear RNA"

		rna_names <- colnames(rna_df)
		rna_names[1] <- "count"
		colnames(rna_df) <- rna_names

		merged <- rbind(circ_df, rna_df)
		merged$count <- log2(merged$count + 1)
		merged_df <- data_summary(merged, varname="count", groupnames=c("condition", "type"))

		make_lineplot(merged_df, inputdata)
		}
}


make_boxplot <- function(circ_df, inputdata){

	circ_id <- inputdata$bed$name

	p1 <- ggboxplot(circ_df, x="condition", y="count",
									fill="condition", palette = "npg",
									title = paste(circ_id),
									ylab = "normalized counts", xlab="",
          				add = c("dotplot"),
									add.params = list(size=0.5, jitter=0.1),
									legend = "none",
									bxp.errorbar = T,
          				bxp.errorbar.width = 0.2, width=0.3,
									ggtheme = theme_classic()) +
          				rotate_x_text(angle = 0) +
									theme(plot.title = element_text(face = "bold", size=16, hjust = 0.5)) +
					  			theme(axis.text.x = element_text( colour = "black", size=14)) +
					  			theme(axis.title.y = element_text(size=14, face = "italic")) +
					  			theme(axis.title.x = element_blank()) +
					  			theme(axis.text.y = element_text(color = "black", size=10))

  pdf(paste(circ_id, "boxplot.pdf", sep="_"))
  plot(p1)
 	dev.off()
}

make_lineplot <- function(merged_df, inputdata){

	circ_id <- inputdata$bed$name
	gene_id <- inputdata$parent

	p2 <- ggline(merged, x="condition", y="count",
				color="type", ylab="log2(normalized counts + 1)",
				title = paste(circ_id, "|", gene_id, sep=" "),
				palette="npg", add = c("mean_se", "jitter"),
				add.params = list(size=0.5, jitter=0.1),
				legend="right", ggtheme = theme_classic()) +
				theme(plot.title = element_text(face = "bold", size=16, hjust = 0.5)) +
  			theme(axis.text.x = element_text( colour = "black", size=14)) +
  			theme(axis.title.y = element_text(size=14, face = "italic")) +
  			theme(axis.title.x = element_blank()) +
  			theme(axis.text.y = element_text(color = "black", size=10))

	pdf(paste(circ_id, gene_id, "expression.pdf", sep="_"))
  plot(p2)
  dev.off()
}




# Packages + Error traceback (really handy for explicit error tracing)
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))

arg <- get_args()
inputdata <- stage_data(arg$de_circ, arg$circ_counts, arg$gene_counts, arg$parent_gene, arg$bed, arg$mature_len, arg$phenotype)
dir.create(inputdata$bed$name)
x <- prep_plots(inputdata)
z <- singular_report(inputdata)
