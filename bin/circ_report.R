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

usage <- function(){giveError("USAGE: circ_report.R de_circ circrna_counts gene_counts bed miranda targetscan mature_len phenotype")}

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

	write.table(out_df, file.path(file_name, paste(file_name, "DESeq2_stats.txt", sep="_")), quote=F, sep="\t", row.names=F)
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
		make_ratioplot(merged_df, inputdata)
		}
}


make_boxplot <- function(circ_df, inputdata){

	circ_id <- inputdata$bed$name

	p1 <- ggplot(circ_df, aes(condition, count, color=condition)) +
  	stat_boxplot(geom = 'errorbar', width=0.1) +
  	geom_boxplot(notch = F, aes(fill = factor(condition)), width=0.3) +
  	theme_bw() +
  	geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +
  	scale_color_manual(values=c('black','black')) +
  	scale_fill_manual(values = c('steelblue1','salmon')) +
  	labs(x = "Condition") +
  	labs(y = "Normalized Counts") +
  	labs(title = paste(circ_id)) +
  	guides(fill = FALSE) +
  	theme(plot.title = element_text(face = "bold", size=18, hjust = 0.5)) +
  	theme(axis.text.x = element_text( colour = "black", size=14)) +
  	theme(axis.title.y = element_text(size=14)) +
  	theme(axis.title.x = element_blank()) +
  	theme(plot.caption = element_text(face = "italic", size=12)) +
  	theme(axis.text.y = element_text(color = "black", size=10))
  	p2 <- p1  + theme(legend.position="none")
  	pdf(paste(circ_id, "Boxplot.pdf", sep="_"))
  	plot(p2)
 	dev.off()
}

make_lineplot <- function(merged_df, inputdata){

	circ_id <- inputdata$bed$name
	gene_id <- inputdata$parent

	p3 <- ggplot(merged_df, aes(x=condition, y=count, group=type, color=type)) +
  	geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=0.1) +
  	geom_line() + geom_point() +
  	scale_color_manual(values=c('dodgerblue1','firebrick1')) + theme_bw() +
  	labs(y = "log2(Normalized Counts)") +
  	labs(title = paste(circ_id, "|", gene_id, sep=" ")) +
  	guides(fill = FALSE) +
  	theme(plot.title = element_text(face = "bold", size=18, hjust = 0.5)) +
  	theme(axis.text.x = element_text( colour = "black", size=14)) +
  	theme(axis.title.y = element_text(size=14)) +
  	theme(axis.title.x = element_blank()) +
  	theme(plot.caption = element_text(face = "italic", size=12)) +
  	theme(axis.text.y = element_text(color = "black", size=10))
  	pdf(paste(circ_id, gene_id, "Expression.pdf", sep="_"))
  	plot(p3)
  	dev.off()
}

make_ratioplot <- function(merged_df, inputdata){

	circ_id <- inputdata$bed$name
	gene_id <- inputdata$parent

	require(dplyr)
  	ratio_df <- merged_df %>%
  	group_by(condition) %>%
  	summarise(ratio = count[type=="circRNA"]/(count[type=="linear RNA"] + count[type=="circRNA"]))
  	detach(package:dplyr)
  	ratio_df$ratio[is.nan(ratio_df$ratio)] <- 0

 	p4 <- ggplot(ratio_df, aes(condition, ratio, fill=condition)) +
  	geom_bar(stat="identity", width=0.4) + # create boxplot
  	theme_bw() +
  	scale_color_manual(values=c('black','black')) +
  	scale_fill_manual(values = c('steelblue1','salmon')) +   # def. col. palette
  	labs(x = "Condition") +
  	labs(y = "Ratio log2(circRNA/circRNA + linear)") +
  	labs(title = paste(circ_id, "|", gene_id, "Ratio", sep=" ")) +
  	guides(fill = FALSE) +
  	theme(plot.title = element_text(face = "bold", size=18, hjust = 0.5)) +
  	theme(axis.text.x = element_text( colour = "black", size=14)) +
  	theme(axis.title.y = element_text(size=14)) +
  	theme(axis.title.x = element_blank()) +
  	theme(plot.caption = element_text(face = "italic", size=12)) +
  	theme(axis.text.y = element_text(color = "black", size=10)) +
  	scale_y_continuous(expand = expansion(mult = c(0, .1)))
  	p5 <- p4  + theme(legend.position="none")
  	pdf(paste(circ_id, "Ratio_Plot.pdf", sep="_"))
  	plot(p5)
  	dev.off()
}



# Packages + Error traceback (really handy for explicit error tracing)
options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("circlize"))

arg <- get_args()
inputdata <- stage_data(arg$de_circ, arg$circ_counts, arg$gene_counts, arg$parent_gene, arg$bed, arg$mature_len, arg$phenotype)
dir.create(inputdata$bed$name)
x <- prep_plots(inputdata)
z <- singular_report(inputdata)
