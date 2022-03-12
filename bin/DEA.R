#!/usr/bin/Rscript

# Automated differential expression analysis script for nf-core/circrna
# Relies on Stirngtie prepde.py outputs
# Would be fun to adapt to STAR + Kallisto for own use in future.

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
            default="phenotype.csv")

    argp <- add_argument(
            parser=argp,
            arg="circRNA",
            short="c",
            help="circRNA counts matrix",
            default="circRNA_matrix.txt")

    argp <- add_argument(
            parser=argp,
            arg="species",
            short="s",
            help="species ID",
            default="hsa")

    argp <- add_argument(
            parser=argp,
            arg="map",
            short="m",
            help="ensDB",
            "default"="ensemblDatabase.txt")

    argv <- parse_args(
            parser=argp,
            argv = commandArgs(trailingOnly = TRUE))

    return(argv)

}

giveError <- function(message){
    cat(paste("\n", message, sep=""))
    quit()
}

usage <- function(){giveError("USAGE: DEA.R <gene_counts.csv> <phenotype.txt> <circRNA_matrix.txt> <species id> <ensembl_map>")}


stage_data <- function(gene_counts, phenotype, circRNA, species, map){

    inputdata <- list()

    dbmap <- read.table(map, sep="\t", header=T, quote="", stringsAsFactors=FALSE)
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
    inputdata$species <- species
    inputdata$map <- dbmap

    inputdata$gene <- ens2symbol(inputdata$gene, inputdata)

    return(inputdata)
}


checkinputdata <- function(phenotype){

    # Stage phenotype file
    pheno <- read.csv(phenotype, row.names=1, header = T, stringsAsFactors = T)

    # Check if there are at least 3 replicates (DESeq2 fails if < 3)
    if(min(table(pheno$condition)) >= 3){
        print("Suitable sample size for DE analysis")
    }else{
        giveError("Not enough samples per condition to perform DE analysis!")
    }

    # Rename sex to gender..
    if("sex" %in% names(pheno)){
        print("Renaming sex to gender in phenotype file")
        rename <- gsub("sex", "gender", names(pheno))
        names(pheno) <- rename
    }

    # Check gender is only male, female, unknown
    if ("gender" %in% names(pheno)) {
        if (! all(unique(pheno$gender) %in% c("m", "f", "u"))) {
            giveError("SAMPLEINFO ERROR:\nOnly the values m [male], f [female] and u [unknown] are supported in field <gender>.\n")
        }
    }

    ## check if all columns are factors. If numeric, convert to factor.
    factor_cols <- sapply(pheno, is.factor)
    if(all(factor_cols) == TRUE){
        print("All columns in phenotype are factors and suitable for analysis.")
    }else{
        numeric_cols <- sapply(pheno, is.numeric)
        names <- colnames(pheno)[numeric_cols]
        print(paste0("Column(s) ", names, " is numeric. Converting to factor."))
        pheno[numeric_cols] <- as.data.frame(lapply(pheno[numeric_cols], factor))
        final_check <- sapply(pheno, is.factor)
        if(all(final_check) == TRUE){
            print("Finished coverting to factor")
        }else{
            giveError("Error in converting to factors. See checkinputdata function.")
        }
    }

    return(pheno)

}



makedesign <- function(phenotype){

    # Covariates i.e explanatory variables.
    covariates <- names(phenotype)[which(!names(phenotype) %in% c("condition"))]
    design <- formula(paste("~", paste(c(covariates, "condition"), sep="", collapse=" + ")))
    return(design)

}




ens2symbol <- function(mat, inputdata){

    ## designed to work on input gene_count_matrix.csv file
    ## everything else downstream no longer needs to be converted

    ## figure out if working with ENS, or ENS IDs

    mat <- as.data.frame(mat)
    map <- inputdata$map
    species <- inputdata$species

    if(all(grepl(pattern="^ENSG", rownames(mat)))){
        filter = "ensembl_gene_id"
            if(all(grepl(pattern=".", rownames(mat)))){
                filter = "ensembl_gene_id_version"
            }
    }else{
        filter = "external_gene_name"
    }

    if(filter == "external_gene_name"){
        print("Using external gene name as gene symbol")
    }else{
        print("Setting up Mart to convert ENS IDs to gene symbols")
        ## set up Mart
        mart_call <- as.character(subset(map$command, map$species == species))
        print("ENS2SYMBOL")
        mart <- eval(str2expression(mart_call))
    }

    ## now go about converting ENS2SYMBOL
    if(filter == "ensembl_gene_id"){

        mat$ensembl_gene_id <- rownames(mat)
        info <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                        filters = c("ensembl_gene_id"),
                        values = mat$ensembl_gene_id,
                        mart = mart,
                        useCache=FALSE)

        tmp <- merge(mat, info, by="ensembl_gene_id")
        tmp$external_gene_name <- make.names(tmp$external_gene_name, unique = T)
        rownames(tmp) <- tmp$external_gene_name
        tmp <- subset(tmp, select=-c(ensembl_gene_id, external_gene_name))

        mat <- tmp
        print("input mat ensembl gene id detected and converted")
        return(mat)
    }else if(filter == "ensembl_gene_id_version"){

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

        mat <- tmp
        print("input mat ensembl gene id version detected and converted")
        return(mat)
    }else{
        print("NO change made to input mat ")
        return(mat)
    }

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


annotate_de_genes <- function(df, inputdata){

    map <- inputdata$map
    species <- inputdata$species
    print("ANNOTATE DE GENES")
    mart_call <- as.character(subset(map$command, map$species == species))
    mart <- eval(str2expression(mart_call))

    df$external_gene_name <- rownames(df)
    info <- getBM(attributes=c("external_gene_name",
                                "chromosome_name",
                                "start_position",
                                "end_position",
                                "strand",
                                "entrezgene_description"),
                    filters = c("external_gene_name"),
                    values = rownames(df),
                    mart = mart,
                    useCache=FALSE)

    tmp <- merge(df, info, by="external_gene_name")
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

# Data type provided at end of script to activate RNA-Seq / circRNA analysis.
DESeq2 <- function(inputdata, data_type){

    if(data_type == "RNA-Seq"){
        outdir <- "RNA-Seq/"

        dds <- DESeqDataSetFromMatrix(
        countData=inputdata$gene,
        colData=inputdata$pheno,
        design = inputdata$design)

        levels <- as.character(unique(inputdata$pheno$condition))
        for(level in levels){
            reference <- level
            contrasts <- levels[levels != paste0(reference)]
            dds$condition <- relevel(dds$condition, ref = paste0(reference))
            dds <- DESeq(dds, quiet=TRUE)

            DESeq2_plots(dds, outdir)

            for(var in contrasts){
                contrast <- paste(var, "vs", reference, sep="_")
                DEG <- getDESeqDEAbyContrast(dds, contrast, reference, var, outdir, inputdata)
            }
        }
    }else if(data_type == "circRNA"){
        outdir <- "circRNA/"

        ## use gene sizeFactors
        tmp <- DESeqDataSetFromMatrix(
        countData=inputdata$gene,
        colData=inputdata$pheno,
        design = inputdata$design)
        tmp <- DESeq(tmp, quiet=TRUE)

        sizefactors <- sizeFactors(tmp)
        rm(tmp)

        dds <- DESeqDataSetFromMatrix(
        countData=inputdata$circ,
        colData=inputdata$pheno,
        design = inputdata$design)

        levels <- as.character(unique(inputdata$pheno$condition))
        for(level in levels){
            reference <- level
            contrasts <- levels[levels != paste0(reference)]
            dds$condition <- relevel(dds$condition, ref = paste0(reference))
            dds <- DESeq(dds, quiet=TRUE)
            sizeFactors(dds) <- sizefactors

            DESeq2_plots(dds, outdir)

            for(var in contrasts){
                contrast <- paste(var, "vs", reference, sep="_")
                DEG <- getDESeqDEAbyContrast(dds, contrast, reference, var, outdir)
            }
        }
    }else{
        giveError("Data type not provided correctly, check end of script")
    }
    return(DEG)
}


getDESeqDEAbyContrast <- function(dds, contrast, reference, var, outdir, inputdata) {

    res <- results(dds, filterFun=ihw, alpha=0.05,  contrast=c("condition", var, reference))
    cat('\n\nSummary data from DESeq2 for ', contrast, ':', sep="")
    summary(res)

    ma_plot(res, contrast, outdir)

    up_regulated <- get_upregulated(res)
    down_regulated <- get_downregulated(res)

    de_up <- rownames(up_regulated)
    de_down <- rownames(down_regulated)
    de <- c(de_up, de_down)
    cts <- counts(dds, normalized=T)

    # attempt boxplots here
    if(outdir == "circRNA/"){
        make_boxplots(de, cts, contrast)
    }

    log2 <- log2(cts +1)
    global_heatmap(de, log2, contrast, outdir)

    if(outdir == "RNA-Seq/"){
        up_regulated <- annotate_de_genes(up_regulated, inputdata)
        down_regulated <- annotate_de_genes(down_regulated, inputdata)
    }else{
        up_regulated <- tibble::rownames_to_column(up_regulated, "ID")
        down_regulated <- tibble::rownames_to_column(down_regulated, "ID")
    }

    dir <- paste(outdir, contrast, sep="")
    dir.create(dir)
    write.table(up_regulated, file.path(dir, paste("DESeq2", contrast, "up_regulated_differential_expression.txt", sep="_")), sep="\t", row.names=F, quote=F)
    write.table(down_regulated, file.path(dir, paste("DESeq2", contrast, "down_regulated_differential_expression.txt", sep="_")), sep="\t", row.names=F, quote=F)

    res_df <- as.data.frame(res)

    #if(outdir == "RNA-Seq/"){
    #    ann_res <- ens2symbol(res_df, inputdata)
    #}else{
    #    ann_res <- res_df
    #}

    volcano_plot(res_df, contrast, outdir)

    pdf(file.path(dir, paste("DESeq2", contrast, "fold_change_distribution.pdf", sep="_")), width=8, height=8)
    hist(res$log2FoldChange, breaks=50, col="seagreen", xlab=paste("(Fold change)", contrast, sep=" "), main="Distribution of differential expression fold change")
    abline(v=c(-1,1), col="black", lwd=2, lty=2)
    legend("topright", "Fold change <-1 and >1", lwd=2, lty=2)
    dev.off()

    pdf(file.path(dir, paste("DESeq2", contrast, "pvalue_distribution.pdf", sep="_")), width=8, height=8)
    hist(res$pvalue, breaks=50, col="seagreen", xlab=paste("P-Value (Fold change)", contrast, sep=" "), main="Distribution of P-Values")
    abline(v=c(0.05),col="black",lwd=2,lty=2)
    legend("topright", "P-Value <0.05",lwd=2,lty=2)
    dev.off()

    pdf(file.path(dir, paste("DESeq2", contrast, "Adj_pvalue_distribution.pdf", sep="_")), width=8, height=8)
    hist(res$padj, breaks=50, col="seagreen", xlab=paste("P-Adj (Fold change)", contrast, sep=" "), main="Distribution of AdjP-Values")
    abline(v=c(0.05),col="black",lwd=2,lty=2)
    legend("top", "P-Adj <0.05",lwd=2,lty=2)
    dev.off()
}


DESeq2_plots <- function(dds, outdir){

    dir.create("DESeq2_QC")
    dir.create(paste("DESeq2_QC/", outdir, sep=""))
    dir=paste("DESeq2_QC/", outdir, sep="")

    pdf(file.path(dir, "DESeq2_dispersion.pdf"), width=8, height=8)
    plotDispEsts(dds)
    dev.off()

    counts <- counts(dds, normalized=T)

    if(outdir == "RNA-Seq/"){
        counts <- ens2symbol(counts, inputdata)
        log2 <- log2(counts + 1)
    }else{
        counts <- as.data.frame(counts)
        log2 <- log2(counts + 1)
    }

    write_counts <- tibble::rownames_to_column(counts, "ID")
    write_log2 <- tibble::rownames_to_column(log2, "ID")

    write.table(write_counts, file.path(outdir, "DESeq2_normalized_counts.txt"), sep="\t", quote=F, row.names = F)
    write.table(write_log2, file.path(outdir, "DESeq2_log2_transformed_counts.txt"), sep="\t", quote=F, row.names = F)

    sample_to_sample_heatmap(log2, outdir)
    sample_to_sample_dendogram(log2, outdir)
    PCA_plot(log2, outdir)

}


ma_plot <- function(res, contrast, outdir){
    dir <- paste(outdir, contrast, sep="")
    dir.create(dir)
    pdf(file.path(dir, paste("DESeq2", contrast, "MA_plot.pdf", sep="_")), width=8, height=8)
    plotMA(res)
    dev.off()

}

sample_to_sample_heatmap <- function(log2, outdir){

    dir.create("DESeq2_QC")
    dir.create(paste("DESeq2_QC/", outdir, sep=""))
    dir=paste("DESeq2_QC/", outdir, sep="")

    sampleDists <- dist(t(log2))
    sampleDistMatrix <- as.matrix(sampleDists)
    pdf(file.path(dir, "DESeq2_sample_heatmap.pdf"), width=8, height=8)
    pheatmap(mat=sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255),
            fontsize_row=8)
    dev.off()

}



sample_to_sample_dendogram <- function(log2, outdir){

    dir.create("DESeq2_QC")
    dir.create(paste("DESeq2_QC/", outdir, sep=""))
    dir=paste("DESeq2_QC/", outdir, sep="")

    d=t(log2)
    d=dist(d)
    hc=hclust(d, method="complete")
    print("test hclust")
    print(head(d))
    pdf(file.path(dir, "DESeq2_sample_dendogram.pdf"))
    plot(hc)
    dev.off()

}


PCA_plot <- function(log2, outdir){

    p <- pca(log2, metadata=inputdata$pheno)

    for(exp_var in names(inputdata$pheno)){
        dir.create("DESeq2_QC")
        dir=paste("DESeq2_QC/", outdir, sep="")
        pdf(file.path(dir, paste("DESeq2", exp_var, "PCA.pdf", sep="_")))
        biplot <- biplot(p,
                        colby=paste(exp_var),
                        hline=0,
                        vline=0,
                        legendPosition="right",
                        legendLabSize=12,
                        legendIconSize=8,
                        lab = TRUE,
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
    symmetric_plot <- max(max_width, abs(min_width))
    min_width <- symmetric_plot * -1
    max_width <- symmetric_plot
    max_height <- -log10(min(res[res$pvalue>0, 5]))

    up <- subset(res, res$log2FoldChange > 1 & res$pvalue <= 0.05)
    up <- up[order(-up$log2FoldChange),]
    up_list <- head(rownames(up), n=10L)

    down <- subset(res, res$log2FoldChange < 1 & res$pvalue <= 0.05)
    down <- down[order(down$log2FoldChange),]
    down_list <- head(rownames(down), n=10L)

    plot_top_20 <- c(up_list, down_list)

    dir <- paste(outdir, contrast, sep="")
    dir.create(dir)
    pdf(file.path(dir, paste("DESeq2", contrast, "volcano_plot.pdf", sep="_")))
    p <- EnhancedVolcano(res,
                        lab=rownames(res),
                        x="log2FoldChange",
                        y="pvalue",
                        selectLab=FALSE,
                        drawConnectors=FALSE,
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


global_heatmap <- function(de, log2, contrast, outdir){

    # Split contrast e.g normal_vs_tumor -> "normal", "tumor"
    levels <- unlist(strsplit(contrast, "_vs_"))
    pheno <- inputdata$pheno
    # subset phenotype file for contrast samples
    pheno_subset <- subset(pheno, pheno$condition %in% levels)
    # check it worked?
    #print(pheno_subset)
    # subset log2 counts for contrast samples
    mat <- log2[,rownames(pheno_subset)]
    # subset de genes/circRNAs
    mat <- mat[de,]
    # Perform scaling and centering on DE expr data
    mat <- t(mat)
    mat <- scale(mat, center=T)
    mat <- t(mat)

    dir <- paste(outdir, contrast, sep="")
    dir.create(dir)
    pdf(file.path(dir, paste("DESeq2", contrast, "heatmap.pdf", sep="_")))
    pheatmap(mat,
            annotation_col=pheno_subset,
            color=greenred(75),
            cluster_rows = T,
            show_rownames = F)
    dev.off()

}



make_boxplots <- function(de, cts, contrast){

    pheno <- inputdata$pheno
    levels <- unlist(strsplit(contrast, "_vs_"))
    pheno_subset <- subset(pheno, pheno$condition %in% levels)
    # subset counts for levels of interest
    counts <- cts[,rownames(pheno_subset)]
    # subset for de genes
    counts <- as.data.frame(counts[de,])
    dir.create("boxplots")
    dir.create(paste("boxplots/", contrast, sep=""))
    dir=paste("boxplots/", contrast, sep="")
    for( i in 1:nrow(counts)){
        circ_id <- rownames(counts[i,]);
        mat <- as.data.frame(t(counts[i,]));
        mat <- cbind(mat, pheno_subset);
        names <- c("counts", "condition");
        colnames(mat) <- names;

        p1 <- ggboxplot(mat, x="condition", y="counts",
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

        pdf(file.path(dir, paste(circ_id, "boxplot.pdf", sep="_")))
        plot(p1)
        dev.off()
    }
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
suppressPackageStartupMessages(library("ggpubr"))
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

inputdata <- stage_data(arg$gene_counts, arg$phenotype, arg$circRNA, arg$species, arg$map)
dir.create("RNA-Seq")
dir.create("circRNA")
x <- DESeq2(inputdata, "RNA-Seq")
y <- DESeq2(inputdata, "circRNA")
