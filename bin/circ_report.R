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
		arg="mature_len",
		short="l",
		help="length of mature circRNA")

	argp <- add_argument(
		parser=argp,
		arg="phenotype",
		short="p",
		help="Phenotype / colData / samples file used previously for DESeq2")

	argp <- add_argument(
		parser=argp,
		arg="circlize_exons",
		short="z",
		help="output from bash script preparing circlize exons")

    	argv <- parse_args(
            	parser=argp,
            	argv = commandArgs(trailingOnly = TRUE))

    return(argv)
    }

giveError <- function(message){
    cat(paste("\n", message, sep=""))
    quit()
    }

usage <- function(){giveError("USAGE: circ_report.R de_circ circrna_counts gene_counts bed miranda targetscan mature_len")}

stage_data <- function(de_circ, circ_counts, gene_counts, parent_gene, bed, miranda, targetscan, mature_len, phenotype, circlize_exons){

	inputdata <- list()

	de <- read.table(de_circ, sep="\t", row.names="ID", header=T, stringsAsFactors=F)
	circ <- read.table(circ_counts, sep="\t", row.names="ID", header=T)
	gene <- read.table(gene_counts, sep="\t", row.names="ID", header=T, check.names=F)
	parent_tmp <- read.table(parent_gene, sep="\t", row.names=1, header=F, stringsAsFactors=F)
	colnames(parent_tmp) <- "gene"
	parent <- parent_tmp$gene
	bed <- read.table(bed, sep="\t", header=F, stringsAsFactors=F)
	colnames(bed) <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "type", "ExonCount", "ExonSizes", "ExonStart")
	miranda <- read.table(miranda, sep="\t", header=T, stringsAsFactors=F)
	targetscan <- read.table(targetscan, sep="\t", header=T, stringsAsFactors=F)
	mature_tmp <- read.table(mature_len)
	mature <- mature_tmp$V1
	pheno <- read.table(phenotype, sep="\t", header=T, row.names=1)
	circlize_exons <- read.table(circlize_exons, sep="\t")

	inputdata$de <- de
	inputdata$circ <- circ
	inputdata$gene <- gene
	inputdata$parent <- parent
	inputdata$bed <- bed
	inputdata$miranda <- miranda
	inputdata$targetscan <- targetscan
	inputdata$mature <- mature
	inputdata$pheno <- pheno
	inputdata$circlize_exons <- circlize_exons

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

	write.table(out_df, file.path(file_name, paste(file_name, "Report.txt", sep="_")), quote=F, sep="\t", row.names=F)
} 
	


miRNAs <- function(inputdata){
	
	outdir <- inputdata$bed$name

	miranda <- inputdata$miranda
	targetscan <- inputdata$targetscan
	x <- inputdata$circlize_exons
	circlize_exons <- inputdata$circlize_exons

	miranda$miRNA <- gsub("hsa-", "", miranda$miRNA)
	colnames(targetscan)[2] <- "miRNA"
	
	mir_df <- merge(miranda, targetscan, by="miRNA")
	
	## Filtering Step:
	## Removes miRNA <= -20.00Kcal/Mol
	## Removes duplicate miRNAs sharing the same bind site (duplicate miRNA IDs)
	## Keeps the miRNA with the higher "Score"
	mir_df <- subset(mir_df, mir_df$Energy_KcalMol <= -20.00)
	mir_df <- mir_df[order(mir_df$MSA_start, -abs(mir_df$Score)),]
	mir_df <- mir_df[!duplicated(mir_df$MSA_start),]
	miRs <- subset(mir_df, select=c(miRNA, MSA_start, MSA_end))
	colnames(miRs) <- c("miRNA", "Start", "End")
	
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
	
	write.table(write_mirs, file.path(outdir, paste(outdir, "miRNA_targets.txt", sep="_")), sep="\t", row.names=F, quote=F)
	make_circos_plot(inputdata, circlize_exons, circlize_mirs)
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

	print(circ_df)
	print(rna_df)
	
	circ_df <- cbind(circ_df, pheno)
	rna_df <- cbind(rna_df, pheno)

	circ_df$type <- "circRNA"
	rna_df$type <- "linear RNA"

	rna_names <- colnames(rna_df)
	rna_names[1] <- "count"
	colnames(rna_df) <- rna_names

	circ_names <- colnames(circ_df)
	circ_names[1] <- "count"
	colnames(circ_df) <- circ_names

	make_boxplot(circ_df, inputdata)

	merged <- rbind(circ_df, rna_df)
	merged$count <- log2(merged$count + 1)
	merged_df <- data_summary(merged, varname="count", groupnames=c("condition", "type"))

	make_lineplot(merged_df, inputdata)
	make_ratioplot(merged_df, inputdata)
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
  	pdf(file.path(circ_id, paste(circ_id, "Boxplot.pdf", sep="_")))
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
  	pdf(file.path(circ_id, paste(circ_id, gene_id, "Expression.pdf", sep="_")))
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
  	pdf(file.path(circ_id, paste(circ_id, "Ratio_Plot.pdf", sep="_")))
  	plot(p5)
  	dev.off()
}

make_circos_plot <- function(inputdata, circlize_exons, circlize_mirs){

	circ_id <- inputdata$bed$name
	col_text <- "grey25"

	pdf(file.path(circ_id, paste(circ_id, "miRNA_Plot.pdf", sep="_")))
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
