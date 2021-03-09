#!/usr/bin/Rscript

##
##
## Pass to script:
## 1. gene_counts matrix from prep_DE.py
## 2. sample file that is essentially colData. Will deduce Design / model matrix from this
## 3. circRNA matrix
## 4. output dir (./) for nextflow
##
##


get_args <- function(){

	argp <- arg_parser(
		description="this script is designed to automate DESeq2 analysis for circRNA nextflow script",
		hide.opts=TRUE)

	argp <- add_argument(
		parser=argp,
		arg="gene_counts",
		short="g",
		help="gene_count_matrix.csv produced by prepDE.py downstream of Stringtie quant with -e flag",
		default="gene_count_matrix.csv")

	argp <- add_argument(
            	parser=argp,
            	arg="phenotype",
		short="p",
            	help="file containing sample metadata information to guide the design",
		default="phenotype.txt")

 	argp <- add_argument(
            	parser=argp,
            	arg="circRNA",
		short="c",
            	help="circRNA counts matrix",
		default="circRNA_counts.txt")

    	argv <- parse_args(
            	parser=argp,
            	argv = commandArgs(trailingOnly = TRUE))

    return(argv)
    }

giveError <- function(message){
    cat(paste("\n", message, sep=""))
    quit()
    }

usage <- function(){giveError("USAGE: DEA.R <gene_counts.csv> <colData.txt> <outDir> ")}


stage_data <- function(gene_counts, phenotype, circRNA){

	inputdata <- list()

	gene_mat <- read.csv(gene_counts, row.names="gene_id", check.names=F)
	circ <- read.table(circRNA, sep ="\t", header = T, stringsAsFactors=FALSE)

	# Merge circRNA genomic loci to ID
	circ$circ <- with(circ, paste0(Chr, sep=":", Start, sep="-", Stop, sep=":", Strand))
	rownames(circ) <- circ$circ
	circ <- subset(circ, select=-c(Chr, Start, Stop, Strand, circ))

	## add pseudocount of 1
	gene_mat <- gene_mat + 1
	circ <- circ + 1

	inputdata$pheno <- checkinputdata(phenotype)

	# make covariates factors (partic if numeric)
	factor_me <- colnames(inputdata$pheno)
	inputdata$pheno[factor_me] <- lapply(inputdata$pheno[factor_me], factor)

	cols <- rownames(inputdata$pheno)

	if(identical(rownames(inputdata$pheno), colnames(gene_mat))){
		circ <- circ[, cols,]
	}else{
		giveError(c("Samples in phenotype file do not match sequencing sample names.\n",
			    "Please check that phenotype samples match gene_count_matrix.csv headers.\n",
			    "*Make sure they are sorted in alphabetical order:",
			    "tail -n +2 phenotype.txt | sort -k1,1\n\n"))
	}

        inputdata$gene <- gene_mat
        inputdata$circ <- circ
	inputdata$design <- makedesign(inputdata$pheno)

	return(inputdata)
}


checkinputdata <- function(phenotype){

	pheno <- read.table(phenotype, sep="\t", row.names=1, header = T)

	if(min(table(pheno$condition)) >= 3){
		print("Suitable sample size for DE analysis")
	}else{
		giveError("Not enough samples per condition to perform DE analysis!")
	}

	if("sex" %in% names(pheno)){
		print("Renaming sex to gender in phenotype file")
		rename <- gsub("sex", "gender", names(pheno))
		names(pheno) <- rename
	}

	if ("gender" %in% names(pheno)) {
        if (! all(unique(pheno$gender) %in% c("m", "f", "u"))) {
        	giveError("SAMPLEINFO ERROR:\nOnly the values m [male], f [female] and u [unknown] are supported in field <gender>.\n")}
        }


	return(pheno)
}



makedesign <- function(phenotype){
	#phenotype = inputdata$pheno
	covariates <- names(phenotype)[which(!names(phenotype) %in% c("condition"))]
	design <- formula(paste("~", paste(c(covariates, "condition"), sep="", collapse=" + ")))
	return(design)
}




ens2symbol <- function(mat){

	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

	mat <- as.data.frame(mat)
	mat$ensembl_gene_id_version <- rownames(mat)

	info <- getBM(attributes=c("ensembl_gene_id_version","external_gene_name"),
	              filters = c("ensembl_gene_id_version"),
		      values = mat$ensembl_gene_id_version,
		      mart = mart,
          useCache=FALSE)

	tmp <- merge(mat, info, by="ensembl_gene_id_version")
	tmp$external_gene_name <- make.names(tmp$external_gene_name, unique = T)
	rownames(tmp) <- tmp$external_gene_name
	tmp <- subset(tmp, select=-c(ensembl_gene_id_version, external_gene_name))

	return(tmp)
}

get_upregulated <- function(df){

	key <- intersect(rownames(df)[which(df$log2FoldChange>=1)], rownames(df)[which(df$pvalue<=0.05)])

  	results <- as.data.frame((df)[which(rownames(df) %in% key),])
	return(results)
}

get_downregulated <- function(df){

  	key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)], rownames(df)[which(df$pvalue<=0.05)])

  	results <- as.data.frame((df)[which(rownames(df) %in% key),])
  	return(results)
}


annotate_de_genes <- function(df){

	df$ensembl_gene_id_version <- rownames(df)

	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

	info <- getBM(attributes=c("ensembl_gene_id_version",
                             "external_gene_name",
                             "chromosome_name",
                             "start_position",
                             "end_position",
                             "strand",
                             "entrezgene_description"),
                filters = c("ensembl_gene_id_version"),
                values = df$ensembl_gene_id_version,
                mart = mart,
                useCache=FALSE)


	tmp <- merge(df, info, by="ensembl_gene_id_version")
	tmp$strand <- gsub("-1", "-", tmp$strand)
	tmp$strand <- gsub("1", "+", tmp$strand)
	tmp$external_gene_name <- make.names(tmp$external_gene_name, unique = T)

	output_col <- c("Gene", "Chromosome", "Start", "Stop", "Strand", "Description", "Log2FC", "P-value", "Adj P-value")
	tmp <- subset(tmp, select=c(external_gene_name, chromosome_name, start_position, end_position, strand, entrezgene_description, log2FoldChange, pvalue, padj))
	colnames(tmp) <- output_col

	if(min(tmp$Log2FC) > 0){
		tmp <- tmp[order(-tmp$Log2FC),]
	}else{
		tmp <- tmp[order(tmp$Log2FC),]
	}

	return(tmp)
}

DESeq2 <- function(inputdata, data_type){

	if(data_type == "RNA-Seq"){

		outdir <- "RNA-Seq/"

		dds <- DESeqDataSetFromMatrix(
		countData=inputdata$gene,
		colData=inputdata$pheno,
		design = inputdata$design)

		dds$condition <- relevel(dds$condition, ref="control")
		dds <- DESeq(dds, quiet=TRUE)
		contrasts <- levels(dds$condition)

	}else if(data_type == "circRNA"){

		outdir <- "circRNA/"

		## use gene sizeFactors

		tmp <- DESeqDataSetFromMatrix(
		countData=inputdata$gene,
		colData=inputdata$pheno,
		design = inputdata$design)
		tmp$condition <- relevel(tmp$condition, ref="control")
		tmp <- DESeq(tmp, quiet=TRUE)
		contrasts <- levels(tmp$condition)

		sizefactors <- sizeFactors(tmp)

		dds <- DESeqDataSetFromMatrix(
		countData=inputdata$circ,
		colData=inputdata$pheno,
		design = inputdata$design)

		dds$condition <- relevel(dds$condition, ref="control")
		dds <- DESeq(dds, quiet=TRUE)
		contrasts <- levels(dds$condition)
		sizeFactors(dds) <- sizefactors

	}else{
		print("ERROR in Execution stage of script")
	}

	DESeq2_plots(dds, outdir)

	for(group in contrasts[contrasts != "control"]){

		DEG <- getDESeqDEAbyContrast(dds, group, outdir)
	}

	return(DEG)
}


getDESeqDEAbyContrast <- function(dds, group, outdir) {
	contrast <- paste("control_vs_", group, sep="")
    	res <- results(dds, filterFun=ihw, alpha=0.05,  contrast=c("condition", group, "control"))
    	cat('\n\nSummary data from DESeq2 for ', contrast, ':', sep="")
    	summary(res)

	up_regulated <- get_upregulated(res)
	down_regulated <- get_downregulated(res)

	de_up <- rownames(up_regulated)
        de_down <- rownames(down_regulated)
        de <- c(de_up, de_down)
        cts <- counts(dds, normalized=T)
        log2 <- log2(cts +1)
        global_heatmap(de, log2, contrast, group, outdir)

	if(outdir == "RNA-Seq/"){
		up_regulated <- annotate_de_genes(up_regulated)
		down_regulated <- annotate_de_genes(down_regulated)
	}else{
		up_regulated <- tibble::rownames_to_column(up_regulated, "ID")
		down_regulated <- tibble::rownames_to_column(down_regulated, "ID")

	}

	write.table(up_regulated, file.path(outdir, paste("DESeq2", contrast, "up_regulated_differential_expression.txt", sep="_")), sep="\t", row.names=F, quote=F)
	write.table(down_regulated, file.path(outdir, paste("DESeq2", contrast, "down_regulated_differential_expression.txt", sep="_")), sep="\t", row.names=F, quote=F)

	res_df <- as.data.frame(res)

	if(outdir == "RNA-Seq/"){
		ann_res <- ens2symbol(res_df)
	}else{
		ann_res <- res_df
	}

	print("starting volcano plot")
	volcano_plot(ann_res, contrast, outdir)

	pdf(file.path(outdir, paste("DESeq2", contrast, "fold_change_distribution.pdf", sep="_")), width=8, height=8)
	hist(res$log2FoldChange, breaks=50, col="seagreen", xlab=paste("(Fold change)", contrast, sep=" "), main="Distribution of differential expression fold change")
	abline(v=c(-1,1), col="black", lwd=2, lty=2)
	legend("topright", "Fold change <-1 and >1", lwd=2, lty=2)
	dev.off()

	pdf(file.path(outdir, paste("DESeq2", contrast, "pvalue_distribution.pdf", sep="_")), width=8, height=8)
	hist(res$pvalue, breaks=50, col="seagreen", xlab=paste("P-Value (Fold change)", contrast, sep=" "), main="Distribution of P-Values")
	abline(v=c(0.05),col="black",lwd=2,lty=2)
	legend("topright", "P-Value <0.05",lwd=2,lty=2)
	dev.off()

	pdf(file.path(outdir, paste("DESeq2", contrast, "Adj_pvalue_distribution.pdf", sep="_")), width=8, height=8)
	hist(res$padj, breaks=50, col="seagreen", xlab=paste("P-Adj (Fold change)", contrast, sep=" "), main="Distribution of AdjP-Values")
	abline(v=c(0.05),col="black",lwd=2,lty=2)
	legend("top", "P-Adj <0.05",lwd=2,lty=2)
	dev.off()
}


DESeq2_plots <- function(dds, outdir){

	pdf(file.path(outdir, "DESeq2_MAplot.pdf"), height=8, width=8)
	plotMA(dds)
	dev.off()

	pdf(file.path(outdir, "DESeq2_dispersion.pdf"), width=8, height=8)
	plotDispEsts(dds)
	dev.off()

	counts <- counts(dds, normalized=T)

	if(outdir == "RNA-Seq/"){
		counts <- ens2symbol(counts)
		log2 <- log2(counts + 1)
	}else{
		counts <- as.data.frame(counts)
		log2 <- log2(counts + 1)
	}

	write_counts <- tibble::rownames_to_column(counts, "ID")
	write_log2 <- tibble::rownames_to_column(log2, "ID")

	write.table(write_counts, file.path(outdir, "DESeq2_normalized_counts.txt"), sep="\t", quote=F, row.names = F)
	write.table(write_log2, file.path(outdir, "DESeq2_log2_transformed_counts.txt"), sep="\t", quote=F, row.names = F)

	print("starting sample heatmaps")
	sample_to_sample_heatmap(log2, outdir)
	print("starting hclust")
	sample_to_sample_dendogram(log2, outdir)
	print("starting PCA")
	PCA_plot(log2, outdir)

}


sample_to_sample_heatmap <- function(log2, outdir){

	sampleDists <- dist(t(log2))
   	sampleDistMatrix <- as.matrix(sampleDists)
    	print("test sample_to_sample_heatmap")
	print(head(sampleDistMatrix))
    	pdf(file.path(outdir, "DESeq2_sample_heatmap.pdf"), width=8, height=8)
    	pheatmap(
	mat=sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_cols=sampleDists,
        col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255),
        fontsize_row=8)
	dev.off()
}



sample_to_sample_dendogram <- function(log2, outdir){

	d=t(log2)
	d=dist(d)
	hc=hclust(d, method="complete")
	print("test hclust")
	print(head(d))
	pdf(file.path(outdir, "DESeq2_sample_dendogram.pdf"), width=8, height=8)
	plot(hc)
	dev.off()
}


PCA_plot <- function(log2, outdir){

	p <- pca(log2, metadata=inputdata$pheno)
  	n_comp <- length(p$components)
  	pdf(file.path(outdir, "DESeq2_Scree_plot.pdf"), width=10, height=8)
  	scree <- screeplot(
		 p,
        	 components = getComponents(p, 1:n_comp),
        	 hline = 80,
		 subtitle="80% variation explained")
	plot(scree)
  	dev.off()

  	for(exp_var in names(inputdata$pheno)){
    		pdf(file.path(outdir, paste("DESeq2", exp_var, "PCA.pdf", sep="_")))
    		biplot <- biplot(
			  p,
       			  colby=paste(exp_var),
			  hline=0,
        		  vline=0,
        		  legendPosition="right",
        		  legendLabSize=12,
        		  legendIconSize=8,
			  lab = FALSE,
			  labSize = 0.0,
        		  drawConnectors=FALSE,
        		  title="PCA bi-plot",
        		  subtitle="PC1 vs. PC2")
		plot(biplot)
    		dev.off()
  	}
}


volcano_plot <- function(res, contrast, outdir){

	res <- na.omit(res)

	min_width <- min(res$log2FoldChange)
	max_width <- max(res$log2FoldChange)
	max_height <- -log10(min(res[res$pvalue>0, 5]))

	up <- subset(res, res$log2FoldChange > 1 & res$pvalue <= 0.05)
	up <- up[order(-up$log2FoldChange),]
	up_list <- head(rownames(up), n=10L)

	down <- subset(res, res$log2FoldChange < 1 & res$pvalue <= 0.05)
	down <- down[order(down$log2FoldChange),]
	down_list <- head(rownames(down), n=10L)

	plot_top_20 <- c(up_list, down_list)
	print("test volcano plot")
	print(min_width)
	print(max_width)
	print(max_height)
	pdf(file.path(outdir, paste("DESeq2", contrast, "volcano_plot.pdf", sep="_")))
	p <- EnhancedVolcano(res,
			x="log2FoldChange",
			y="pvalue",
			FCcutoff=1.0,
			pCutoff=0.05,
			title="Volcano Plot",
			subtitle=paste(contrast),
			legendVisible=F,
			caption = paste0('Total Genes = ', nrow(res)),
			xlim=c(min_width, max_width),
			ylim=c(0, max_height),
			pointSize = 1.5)
	plot(p)
	dev.off()
}


global_heatmap <- function(de, log2, contrast, group, outdir){

				levels <- unlist(strsplit(contrast, "_vs_"))
				pheno <- inputdata$pheno
				pheno_subset <- subset(pheno, pheno$condition %in% levels)
				mat <- log2[,rownames(pheno)]
				mat <- mat[de,]
        #pheno <- inputdata$pheno
        #pheno_mtx <- as.matrix(pheno)
        #index <- which(matrix(grepl(paste(group), pheno_mtx), ncol=ncol(pheno_mtx)), arr.ind=T)
        #index_col <- unique(index[,2])
        #col_anno <- names(pheno[index_col])
        #sample_col <- pheno[paste(col_anno)]
        #mat <- as.data.frame((log2)[which(rownames(log2) %in% de),])
        mat <- t(mat)
        mat <- scale(mat, center=T)
        mat <- t(mat)
        pdf(file.path(outdir, paste("DESeq2", contrast, "heatmap.pdf", sep="_")))
        pheatmap(mat,
                annotation_col=pheno_subset,
                color=greenred(75),
                cluster_rows = T,
                show_rownames = F)
        dev.off()
}


options(error=function()traceback(2))
suppressPackageStartupMessages(library("argparser"))
#suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("dplyr"))
#suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("EnhancedVolcano"))
#suppressPackageStartupMessages(library("EnsDb.Hsapiens.v86"))
#suppressPackageStartupMessages(library("genefilter")) #for rowVars
suppressPackageStartupMessages(library("ggplot2"))
#suppressPackageStartupMessages(library("ggrepel"))
#suppressPackageStartupMessages(library("ggfortify"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("IHW"))
#suppressPackageStartupMessages(library("limma"))
#suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
#suppressPackageStartupMessages(library("readr"))
#suppressPackageStartupMessages(library("Rsubread"))
#suppressPackageStartupMessages(library("tximport"))
#suppressPackageStartupMessages(library("VennDiagram"))

arg <- get_args()

inputdata <- stage_data(arg$gene_counts, arg$phenotype, arg$circRNA)
dir.create("RNA-Seq")
dir.create("circRNA")
x <- DESeq2(inputdata, "RNA-Seq")
y <- DESeq2(inputdata, "circRNA")
