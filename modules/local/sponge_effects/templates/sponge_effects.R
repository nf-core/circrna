#!/usr/bin/env Rscript

library("SPONGE")
library("rtracklayer")
library("doParallel")
library("foreach")
library("dplyr")
library("visNetwork")
library("pheatmap")
library("ggplot2")
library("reshape2")
library("ggpubr")
library("GSVA")

########################
## GENERAL FUNCTIONS ###
########################
# plot ceRNA network with gene annotations
plot_network <- function(ceRNA_network, signif_hits = NULL,
                         gtf = NULL, annotation = NULL) {
  # plot network 
  ceRNA_plot <- sponge_plot_network(ceRNA_network, genes_miRNA_candidates, ) %>%
    visNetwork::visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 1)))
   
  ceRNA_plot\$x\$edges\$label <- paste("mscor:", round(ceRNA_network\$mscor, 2))

  # extract nodes and edges for customization
  nodes <- ceRNA_plot\$x\$nodes
  edges <- ceRNA_plot\$x\$edges
  
  # mark differentially expressed RNAs
  if (!is.null(signif_hits)){
    nodes[nodes\$id %in% hgncs\$hgnc_symbol | nodes\$id %in% signif_hits\$X,"group"] <- "DE"
  }
  if (!is.null(gtf)) {
    gtf <- rtracklayer::readGFF(gtf)
    annotation <- unique(gtf[!is.na(gtf\$transcript_id),c("gene_id", "gene_name", "gene_type")])
    colnames(annotation) <- c("ensembl_gene_id", "hgnc_symbol", "gene_type")
    rownames(annotation) <- annotation\$ensembl_gene_id
  }
  if (!is.null(annotation)){
    # break types of RNA down to lncRNA, coding, and rest
    annotation\$gene_type[annotation\$gene_type != "protein_coding" & annotation\$gene_type != "lncRNA"] <- "other_RNA"
    
    nodes <- merge(nodes, annotation, by = 1, all.x = TRUE)
    
    # change to hgnc
    nodes[!is.na(nodes\$hgnc_symbol), 1:2] <- nodes[!is.na(nodes\$hgnc_symbol), "hgnc_symbol"]
    
    # add circRNA as gene_type
    nodes[is.na(nodes\$gene_type) & !grepl("ENSG", nodes\$id), "gene_type"] <- "circRNA"
    # label unknown RNAs as other
    nodes[is.na(nodes\$gene_type), "gene_type"] <- "other_RNA"
    # remove preset color and shape
    nodes <- nodes[, -c(3,4)]
    # change to group
    colnames(nodes)[8] <- "group"
    
    # convert geneA and geneB
    edges <- merge(edges, annotation, by.x = "from", by.y = 1, all.x = TRUE)
    edges[!is.na(edges\$hgnc_symbol), "from"] <- edges\$hgnc_symbol[!is.na(edges\$hgnc_symbol)]
    edges <- merge(edges, annotation, by.x = "to", by.y = 1, all.x = TRUE)
    edges[!is.na(edges\$hgnc_symbol.y), "to"] <- edges\$hgnc_symbol.y[!is.na(edges\$hgnc_symbol.y)]
  }
  # plot final graph
  graph <- visNetwork(nodes = nodes, edges = edges) %>%
    visIgraphLayout(type = "full", physics = FALSE) %>%
    visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
    visGroups(groupname = "protein_coding", shape = "triangle", color = "#0066CC") %>%
    visGroups(groupname = "lncRNA", shape = "square", color = "#F8766D") %>%
    visGroups(groupname = "other_RNA", color = "#DC71FA", shape = "diamond") %>%
    visGroups(groupname = "DE", color = "#CC3333") %>%
    visLegend()
  return(graph)
}



circ.mRNA.subnetwork <- function(interactions, pattern){
  return(interactions[
    (grepl(pattern, interactions\$geneA) & !grepl(pattern, interactions\$geneB)) |
    (grepl(pattern, interactions\$geneB) & !grepl(pattern, interactions\$geneA)),
  ])
}

# plot model performance (spongEffects)
plot_performance <- function (trained_model, central_genes_model = NULL,
                              random_model, training_dataset_name = "TCGA",
                              testing_dataset_name = "TCGA", subtypes) {
  trained.model <- trained_model
  CentralGenes.model <- central_genes_model
  Random.model <- random_model
  training_string <- paste0(training_dataset_name, " (Training)")
  testing_string <- paste0(testing_dataset_name, " (Testing)")
  set_names <- c(training_string, testing_string)
  SpongingActiivty.model <- trained_model

  if (!is.null(central_genes_model)) {
    Accuracy.df <- data.frame(
      Run = rep(set_names, 3),
      Model = c(rep("Modules", 2), rep("Central Genes", 2), rep("Random", 2)),
      Accuracy = c(
        SpongingActiivty.model\$ConfusionMatrix_training[["overall"]][["Accuracy"]],
        SpongingActiivty.model\$ConfusionMatrix_testing[["overall"]][["Accuracy"]],
        CentralGenes.model\$ConfusionMatrix_training[["overall"]][["Accuracy"]],
        CentralGenes.model\$ConfusionMatrix_testing[["overall"]][["Accuracy"]],
        Random.model\$ConfusionMatrix_training[["overall"]][["Accuracy"]],
        Random.model\$ConfusionMatrix_testing[["overall"]][["Accuracy"]]
      )
    )
    Accuracy.df\$Model <- factor(Accuracy.df\$Model, levels = c("Modules", "Random", "Central Genes"))
    Accuracy.df\$Run <- factor(Accuracy.df\$Run, levels = set_names)
  }
  else {
    Accuracy.df <- data.frame(
      Run = rep(set_names, 2),
      Model = c(rep("Modules", 2),
      rep("Random", 2)),
      Accuracy = c(
        SpongingActiivty.model\$ConfusionMatrix_training[["overall"]][["Accuracy"]],
        SpongingActiivty.model\$ConfusionMatrix_testing[["overall"]][["Accuracy"]],
        Random.model\$ConfusionMatrix_training[["overall"]][["Accuracy"]],
        Random.model\$ConfusionMatrix_testing[["overall"]][["Accuracy"]]
      )
    )
    Accuracy.df\$Model <- factor(Accuracy.df\$Model, levels = c("Modules", "Random"))
    Accuracy.df\$Run <- factor(Accuracy.df\$Run, levels = set_names)
  }
  Accuracy.plot <- Accuracy.df %>% ggplot(aes(x = Accuracy, y = Model))
    + geom_line(aes(group = Model))
    + geom_point(aes(shape = Run))
    + theme_light()
    + xlab("Subset Accuracy")
    + ylab("")
    + theme_bw()
  Metrics.SpongeModules.training <- SpongingActiivty.model\$ConfusionMatrix_training\$byClass[,"Balanced Accuracy"] %>% 
    as.data.frame() %>% mutate(Model = "Modules") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.SpongeModules.training) = c("Class", "Value", "Model")
  Metrics.Random.training <- Random.model\$ConfusionMatrix_training\$byClass[, "Balanced Accuracy"] %>% 
    as.data.frame() %>% mutate(Model = "Random") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.Random.training) = c("Class", "Value", "Model")
  if (!is.null(central_genes_model)) {
    Metrics.CentralGenes.training <- CentralGenes.model\$ConfusionMatrix_training[["byClass"]][c(1:length(unique(subtypes)))] %>%
      as.data.frame() %>% mutate(Model = "Central Genes") %>%
      tibble::rownames_to_column("Class")
    colnames(Metrics.CentralGenes.training) = c("Class", "Value", "Model")
    Metrics.training <- rbind( Metrics.SpongeModules.training, rbind( Metrics.Random.training, Metrics.CentralGenes.training)) %>% mutate(Run = training_string)
  }
  else {
    Metrics.training <- rbind(Metrics.SpongeModules.training, Metrics.Random.training) %>% mutate(Run = training_string)
  }
  Metrics.SpongeModules.testing <- SpongingActiivty.model\$ConfusionMatrix_testing\$byClass[, "Balanced Accuracy"] %>%
    as.data.frame() %>% mutate(Model = "Modules") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.SpongeModules.testing) = c("Class", "Value", "Model")
  Metrics.Random.testing <- Random.model\$ConfusionMatrix_testing\$byClass[, "Balanced Accuracy"] %>%
    as.data.frame() %>% mutate(Model = "Random") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.Random.testing) = c("Class", "Value", "Model")
  if (!is.null(central_genes_model)) {
    Metrics.CentralGenes.testing <- CentralGenes.model\$ConfusionMatrix_testing[["byClass"]][c(1:length(unique(subtypes)))] %>% 
      as.data.frame() %>% mutate(Model = "Central Genes") %>% tibble::rownames_to_column("Class")
    colnames(Metrics.CentralGenes.testing) = c("Class", "Value", "Model")
    Metrics.testing <- rbind(Metrics.SpongeModules.testing, rbind(Metrics.Random.testing, Metrics.CentralGenes.testing)) %>%
      mutate(Run = testing_string)
  }
  else {
    Metrics.testing <- rbind(Metrics.SpongeModules.testing, Metrics.Random.testing) %>% mutate(Run = testing_string)
  }
  Metrics <- rbind(Metrics.training, Metrics.testing)
  if (!is.null(central_genes_model)) {
    Metrics\$Model <- factor(Metrics\$Model, levels = c("Modules", "Random", "Central Genes"))
  }
  else {
    Metrics\$Model <- factor(Metrics\$Model, levels = c("Modules", "Random"))
  }
  Metrics\$Run <- factor(Metrics\$Run, levels = set_names)
  Metrics\$Class <- gsub("Class: ", "", Metrics\$Class)
  Metrics.plot <- Metrics %>% ggplot(aes(x = Class, y = Value, fill = Model))
    + geom_bar(position = "dodge", stat = "identity", width = 0.5)
    + facet_grid(Metrics\$Run)
    + xlab("Accuracy")
    + ylab("") + theme_bw()
  metric_plots <- ggarrange(Accuracy.plot, Metrics.plot, ncol = 1, nrow = 2)
  return(metric_plots)
}

# plot lollipop (spongEffects)
plot_modules <- function (trained_model,
                          k_modules = 25,
                          k_modules_red = 10,
                          text_size = 16) {
  final.model <- trained_model\$Model\$finalModel
  Variable.importance <- importance(final.model) %>% as.data.frame() %>%
    tibble::rownames_to_column("Module") %>% arrange(desc(MeanDecreaseGini))
  grey_modules = k_modules - k_modules_red
  p <- Variable.importance[1:k_modules, ] %>%
    mutate(Analysed = c(rep("1", k_modules_red), rep("0", grey_modules))) %>%
    arrange(desc(MeanDecreaseGini)) %>%
    ggplot(aes(x = reorder(Module, MeanDecreaseGini), y = MeanDecreaseGini))
      + geom_point()
      + geom_segment(aes(x = Module, xend = Module, y = 0, yend = MeanDecreaseGini, color = Analysed))
      + scale_colour_manual(values = c("red", "black"), breaks = c("1", "0"))
      + coord_flip()
      + xlab("Module")
      + ylab("Mean decrease in Gini index")
      + theme_light()
      + theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(), axis.ticks.y = element_blank(),
              legend.title = element_blank(), legend.position = "none",
              legend.background = element_blank(), legend.direction = "horizontal",
              panel.border = element_rect(colour = "black", fill = NA, size = 1), text = element_text(size = 16))
  return(p)
}

# Training heat map (spongEffects)
plot_hmap <- function (trained_model,
                       spongEffects,
                       meta_data,
                       label,
                       sampleIDs,
                       Modules_to_Plot = 5,
                       show.rownames = F,
                       show.colnames = F) {
  if (label %in% colnames(meta_data) & sampleIDs %in% colnames(meta_data)) {
    final.model <- trained_model\$Model\$finalModel
    Variable.importance <- importance(final.model) %>% as.data.frame() %>%
      tibble::rownames_to_column("Module") %>% arrange(desc(MeanDecreaseGini))
    Variable.importance\$Module <- gsub("`", "", Variable.importance\$Module)
    Annotation.meta <- meta_data[match(colnames(spongEffects), meta_data[, sampleIDs]), ]
    unique_subtypes <- unique(Annotation.meta\$label)
    number_groups <- length(unique_subtypes)
    col.heatmap <- met.brewer("Renoir", n = number_groups, type = "continuous")
    col.heatmap <- as.vector(col.heatmap)
    col.heatmap <- setNames(col.heatmap, unique_subtypes)
    col.heatmap <- cell.colors
    Column.Annotation <- HeatmapAnnotation(Group = Annotation.meta[, label], col = list(Group = col.heatmap))
    spongeEffects.toPlot <- spongEffects[match(Variable.importance\$Module, rownames(spongEffects)), ]
    spongeEffects.toPlot <- spongeEffects.toPlot[rowSums(is.na(spongeEffects.toPlot)) == 0, ]
    if (Modules_to_Plot > length(rownames(spongeEffects.toPlot)))
      Modules_to_Plot = length(rownames(spongeEffects.toPlot))
    spongeEffects.toPlot <- spongeEffects.toPlot[1:Modules_to_Plot, ]
    Heatmap.p <- spongeEffects.toPlot %>% t() %>% scale() %>%
      t() %>% Heatmap(show_row_names = show.rownames,
                      show_column_names = show.colnames,
                      top_annotation = Column.Annotation,
                      show_heatmap_legend = TRUE)
    return(Heatmap.p)
  }
  else {
    print("label and/or sampleIDs must be columns in metadata")
  }
}

# Target gene heat map (spongEffects)
plot_target_gene_expressions <- function(target, target_genes,
                                         gene_expression, meta,
                                         log_transform = TRUE,
                                         pseudocount = 1, gtf_raw = NULL,
                                         annotation = NULL,
                                         split = "condition", unit = "counts",
                                         show_rows = TRUE,
                                         annotation_colors = NULL) {
  # get target expression
  target_expression <- gene_expression[ ,target, drop = FALSE]
  # get target genes expressions
  target_genes_expression <- gene_expression[ ,target_genes]
  # combine into one
  data <- t(cbind(target_expression, target_genes_expression))
  # save all conditions of samples
  conditions <- unique(meta[,split])
  
  # convert to hgnc if gtf is given
  if (!is.null(gtf_raw)) {
    print("reading GTF file and converting to HGNC symbols...")
    gtf <- rtracklayer::readGFF(gtf)
    gene.ens.all <- unique(gtf[!is.na(gtf\$transcript_id), c("gene_id", "gene_name")])
    colnames(gene.ens.all) <- c("ensembl_gene_id", "hgnc_symbol")
    rownames(gene.ens.all) <- gene.ens.all\$ensembl_gene_id
    # convert
    annotation = gene.ens.all
  }
  if (!is.null(annotation)){
    ensgs <- rownames(data)
    ensgs <- merge(ensgs, annotation, by = 1, all.x = TRUE)
    ensgs[!is.na(ensgs\$hgnc_symbol),"x"] <- ensgs[!is.na(ensgs\$hgnc_symbol), "hgnc_symbol"]
    rownames(data) <- ensgs\$x
  }
  # heat color scheme
  colors <- c(colorRampPalette(c("blue", "orange"))(100), colorRampPalette(c("orange", "red"))(100))
  # annotation colors are given
  if (!is.null(annotation_colors)) {
    annotation.colors <- annotation_colors
  } else {
    annotation.colors <- hcl.colors(length(conditions), palette = hcl.pals(type = "diverging")[12])
  }
  # name colors
  names(annotation.colors) <- conditions
  a_c <- list(x = annotation.colors)
  names(a_c) <- split

  if (log_transform) {
    data = log2(data + pseudocount)
  }
  colnames(data) <- gsub("\\\\.", "-", colnames(data))
  label <- paste0("Module ", unit, " expression for: ", target)
  # create annotation for heat map
  df <- data.frame(meta[,c("sample", split)], row.names = 1)
  pheatmap::pheatmap(data, treeheight_row = 0, treeheight_col = 0,
                     show_colnames = FALSE, show_rownames = show_rows,
                     cluster_rows = TRUE, cluster_cols = TRUE,
                     color = colors, annotation_col = df,
                     annotation_colors = a_c,
                     main = label, fontsize_row = 10)
}

# plot expression of miRNAs in module for specific RNA (no negative values!)
plot_miRNAs_per_gene <- function(target, genes_miRNA_candidates,
                                 mir_expr, meta,
                                 log_transform = TRUE, pseudocount = 1e-3,
                                 unit = "counts", split = "condition",
                                 annotation_colors = NULL){
  # extract miRNAs associated with given RNA
  candidates <- genes_miRNA_candidates[[target]]
  miRNAs <- candidates[["mirna"]]
  miRNA_expr <- mir_expr[ ,miRNAs]
  
  # sum all expressions for the given conditions
  miRNA_expr <- merge(miRNA_expr, meta[ ,c("sample", split)], by.x = 0, by.y = "sample")
  miRNA_expr <- data.frame(miRNA_expr, row.names = 1, check.names = FALSE)
  miRNA_expr <- aggregate(miRNA_expr[ ,1:(ncol(miRNA_expr)-1)], list(condition = miRNA_expr[ ,split]), FUN=sum)
  miRNA_expr <- t(data.frame(miRNA_expr, row.names = 1))
  miRNA_expr <- melt(miRNA_expr, id.vars = "x", varnames = c("miRNA", "variable"))
  miRNA_expr\$miRNA <- gsub("\\\\.", "-", miRNA_expr\$miRNA)
  
  # label
  label <- paste0("miRNA ", unit, " over samples")
  # log transform if given
  if (log_transform){
    miRNA_expr\$value <- log10(miRNA_expr\$value + pseudocount)
    label <- paste0("log10 + ", pseudocount, " ", label)
  }
  # annotation colors are given
  if (!is.null(annotation_colors)) {
    annotation.colors <- annotation_colors
  } else {
    annotation.colors <- hcl.colors(length(conditions), palette = hcl.pals(type = "diverging")[12])
  }
  # name colors
  names(annotation.colors) <- conditions
  # two bar plots in one
  bars <- ggplot(miRNA_expr, aes(y=miRNA, x=value, fill=variable))
    + ggtitle(target)
    + geom_bar(stat='identity', position='dodge')
    + xlab(label)
    + scale_fill_manual(values = annotation.colors, name = "Condition")
  return(bars)
}

#################################################################
## SPONGE FUNCTIONS DUE TO VERSION CONFLICTS WITH RTRACKLAYSER ##
#################################################################

fn_filter_network <- function(network, mscor.threshold = .1, padj.threshold = .01) {
   network %>% filter(mscor > mscor.threshold & p.adj < padj.threshold)
}

fn_weighted_degree <- function(network, undirected = T, Alpha = 1) {
  # Format input matrix by using numeric as node IDs
  Nodes <- data.frame(Nodes = union(network\$geneA, network\$geneB),
                      Nodes_numeric = seq(1, length(union(network\$geneA, network\$geneB))))

  geneA.numeric <- Nodes\$Nodes_numeric[match(network\$geneA, Nodes\$Nodes)]
  geneB.numeric <- Nodes\$Nodes_numeric[match(network\$geneB, Nodes\$Nodes)]

  Input.network <- data.frame(Sender = geneA.numeric, Receiver = geneB.numeric, Weight = network\$mscor)

  if (undirected) {
    # Define networks as undirected
    Undirected.net <- Input.network %>% tnet::symmetrise_w()
    Weighted_degree <- tnet::degree_w(Undirected.net, alpha = Alpha) %>% as.data.frame()
    Nodes\$Weighted_degree <- Weighted_degree\$output[match(Weighted_degree\$node, Nodes\$Nodes_numeric)]
    return(Nodes)
  }
}

filter_ceRNA_network <- function(sponge_effects, 
                                 Node_Centrality = NA, 
                                 add_weighted_centrality = T, 
                                 mscor.threshold = NA, 
                                 padj.threshold = NA) {

  #Filter SPONGE network for significant edges
  Sponge.filtered <- sponge_effects %>%
    fn_filter_network(mscor.threshold =  mscor.threshold, padj.threshold = padj.threshold)
        
  # Some gene columns from SPONGEdb start with an uppercase G, so we need to adjust it just in case
  if("GeneA" %in% colnames(Sponge.filtered)) {
    Sponge.filtered = Sponge.filtered %>% rename(geneA = GeneA)
  }
  if("GeneB" %in% colnames(Sponge.filtered)) {
    Sponge.filtered = Sponge.filtered %>% rename(geneB = GeneB)
  }
    Node_Centrality <- Node_Centrality %>% dplyr::filter(gene %in% Sponge.filtered\$geneA | gene %in% Sponge.filtered\$geneB)

  # Calculate weighted centrality scores and add them to the ones present in SpongeDB
  if(!add_weighted_centrality) {
    sponge_network_centralites <- list(Sponge.filtered)
    names(sponge_network_centralites) <- c("Sponge.filtered")
  }
  else {
    Nodes <- fn_weighted_degree(Sponge.filtered, undirected = T, Alpha = 1)
    Node_Centrality <- Node_Centrality %>%
       mutate(Weighted_Degree = Nodes\$Weighted_degree[match(Node_Centrality\$gene, Nodes\$Nodes)])
    sponge_network_centralites <- list(Sponge.filtered,Node_Centrality)
    names(sponge_network_centralites) <- c("Sponge.filtered","Node_Centrality")
  }
  return(sponge_network_centralites)
}

#---------------------------OUTPUT DIR--------------------------------
ROOT = "."
# create plot directory
PLOT_DIR = file.path(ROOT, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
#---------------------------SPONGE DATA-------------------------------
print("loading SPONGE data...")
load("${sponge_data}") # input from pipeline

#---------------------------GTF DATA----------------------------------
gtf <- rtracklayer::readGFF("${params.gtf}")
gene.ens.all <- unique(gtf[!is.na(gtf\$transcript_id),c("gene_id", "transcript_id", "gene_name")])
print("loaded gtf")
colnames(gene.ens.all) <- c("ensembl_gene_id", "hgnc_symbol", "gene_type")
rownames(gene.ens.all) <- gene.ens.all\$hgnc_symbol

#---------------------------PARAMETERS--------------------------------
mscor.threshold = 0.001 # ${params.sponge_ef_mscor}
padj.threshold = 0.7 # ${params.sponge_ef_fdr}
modules_cutoff = ${params.sponge_ef_enrichment_modules}
bin.size = ${params.sponge_ef_enrichment_bins}
min.size = ${params.sponge_ef_enrichment_min}
max.size = ${params.sponge_ef_enrichment_max}
min.expr = ${params.sponge_ef_enrichment_expr}
method = "${params.sponge_ef_enrichment_method}"
split_training = ${params.sponge_ef_training}
folds = ${params.sponge_ef_folds}

#---------------------------PARALLEL BACKGROUND-----------------------
num.of.cores <- 8 
cl <- makeCluster(num.of.cores)
registerDoParallel(cl)

#---------------------------META FILE/ CONDITIONS---------------------
meta <- read.csv(file = "${params.input}")

#---------------------------EXTENDING META----------------------------
phenotype <- read.csv("${params.phenotype}")
meta <- merge(meta, phenotype, by = "sample")


# ----------------EXPRESSION SPLITTING--------------------------------
n.train <- round(nrow(meta)*split_training)
n.test <- round(nrow(meta)-n.train)

# randomize samples
meta <- meta[sample(1:nrow(meta)), ]
# take n samples from each group
cond.split <- split(meta, meta\$condition)
cond.ratio <- round(sapply(cond.split, function(x) nrow(x) * split_training))
train.meta <- data.table::rbindlist(mapply(function(x,y) x[1:y,], cond.split, cond.ratio, SIMPLIFY = FALSE))
colnames(train.meta) <- c("sampleIDs", c(colnames(train.meta)[3:length(train.meta) - 1], "label"))
test.meta <- data.table::rbindlist(mapply(function(x,y) x[(y + 1):nrow(x), ], cond.split, cond.ratio, SIMPLIFY = FALSE))
colnames(test.meta) <- c("sampleIDs", c(colnames(test.meta)[3:length(test.meta)], "label"))

rownames(gene_expr) <- gsub("\\\\.", "-", rownames(gene_expr)) # MALTE TODO: check if this translates to \\ in work dir
# train gene expression; split_training of samples
train_gene_expr <- gene_expr[rownames(gene_expr) %in% train.meta\$sampleIDs, ]
train_gene_expr <- t(train_gene_expr)
# test gene expression; rest of samples
test_gene_expr <- gene_expr[rownames(gene_expr) %in% test.meta\$sampleIDs, ]
test_gene_expr <- t(test_gene_expr)

# train miRNA expression
train_mirna_expr <- mir_expr[rownames(mir_expr) %in% train.meta\$sampleIDs, ]
# test miRNA expression
test_mirna_expr <- mir_expr[rownames(mir_expr) %in% test.meta\$sampleIDs, ]

#-------------------------CE RNA SPLITTING------------------------
# MALTE TODO: no circ interactions -> empty df
# MALTE TODO: hard code path for now in template
ceRNA_interactions_sign <- read.table("/nfs/data3/CIRCEST/runs/sponging/sponge_res/ceRNA_interactions_sign.tsv", header = TRUE, check.names = FALSE, sep = "\\t")

ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign\$p.adj < padj.threshold), ]
ceRNA_interactions_fdr <- circ.mRNA.subnetwork(ceRNA_interactions_fdr, "circ")

n.train.ce <- round(nrow(ceRNA_interactions_fdr) * split_training)
n.test.ce <- nrow(ceRNA_interactions_fdr) - n.train.ce
train_ceRNA_interactions <- head(ceRNA_interactions_fdr, n.train.ce)
test_ceRNA_interactions <- tail(ceRNA_interactions_fdr, n.test.ce)
#-------------------------CENTRALITIES SPLITTING------------------
network_centralities  <- sponge_node_centralities(ceRNA_interactions_fdr)
n.train.c <- round(nrow(network_centralities) * split_training)
n.test.c <- nrow(network_centralities) - n.train.c
train_network_centralities <- head(network_centralities, n.train.c)
test_network_centralities <- tail(network_centralities, n.test.c)

#-------------------------FILTERING-------------------------------
print("filtering centralities...")

filtered_network_centralities <- filter_ceRNA_network(sponge_effects = train_ceRNA_interactions,
                                                      Node_Centrality = train_network_centralities,
                                                      add_weighted_centrality = TRUE,
                                                      mscor.threshold = mscor.threshold,
                                                      padj.threshold = padj.threshold)
#-------------------------RNAS_OF_INTEREST------------------------
RNAs <- c("lncRNA", "circRNA")
# ensembl.df is part of SPONGE lib
RNAs.ofInterest <- ensembl.df %>% dplyr::filter(gene_biotype %in% RNAs) %>%
  dplyr::select(ensembl_gene_id)
# add circRNAs of the data set
new_circRNAs <- data.frame(ensembl_gene_id = colnames(gene_expr)[grep("c", colnames(gene_expr))])
RNAs.ofInterest <- rbind(RNAs.ofInterest, new_circRNAs)
# get central modules
central_gene_modules <- get_central_modules(central_nodes = RNAs.ofInterest\$ensembl_gene_id,
                                            node_centrality = filtered_network_centralities\$Node_Centrality,
                                            ceRNA_class = RNAs,
                                            centrality_measure = "Weighted_Degree",
                                            cutoff = modules_cutoff)
#-------------------------SPONGE MODULES---------------------------
Sponge.modules <- define_modules(network = filtered_network_centralities\$Sponge.filtered,
                                 central.modules = central_gene_modules,
                                 remove.central = FALSE,
                                 set.parallel = FALSE)
# Module size distribution
Size.modules <- sapply(Sponge.modules, length)
#-------------------------SPLIT MODULES----------------------------
# train and test modules
train.modules <- enrichment_modules(Expr.matrix = train_gene_expr,
                                    modules = Sponge.modules,
                                    bin.size = bin.size,
                                    min.size = min.size,
                                    max.size = max.size,
                                    min.expr = min.expr,
                                    method = method)

test.modules <-  enrichment_modules(Expr.matrix = test_gene_expr,
                                    modules = Sponge.modules,
                                    bin.size = bin.size,
                                    min.size = min.size,
                                    max.size = max.size,
                                    min.expr = min.expr,
                                    method = method)

#-------------------------MODEL PERFORMANCE--------------------------
common_modules = intersect(rownames(train.modules), rownames(test.modules))
train.modules = train.modules[common_modules, ]
test.modules = test.modules[common_modules, ]
trained.model = calibrate_model(Input = train.modules,
                                modules_metadata = train.meta,
                                label = "label",
                                sampleIDs = "sampleIDs",
                                Metric = "Exact_match",
                                n_folds = folds,
                                repetitions = 3)
trained.model[["ConfusionMatrix_training"]]

Input.test <- t(test.modules) %>% scale(center = TRUE, scale = TRUE)
Prediction.model <- predict(trained.model\$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test.meta\$label))
trained.model\$ConfusionMatrix_testing <- ConfusionMatrix_testing

#-------------------------RANDOM MODULES----------------------------
Random.modules <- Random_spongEffects(sponge_modules = Sponge.modules,
                                      gene_expr = train_gene_expr,
                                      min.size = min.size,
                                      bin.size = bin.size,
                                      max.size = max.size,
                                      min.expression=min.expr,
                                      replace = FALSE,
                                      method = method)

# We can now use the randomly defined modules to calculate their enrichment in the test set
Random.modules.test <- enrichment_modules(Expr.matrix = test_gene_expr,
                                          modules = Random.modules\$Random_Modules,
                                          bin.size = bin.size,
                                          min.size = min.size,
                                          max.size = max.size,
                                          min.expr = min.expr,
                                          method = method)

# We find random modules that were identified both in the train and test and use those as input features for the model
common_modules_random = intersect(rownames(Random.modules\$Enrichment_Random_Modules), rownames(Random.modules.test))
Random.modules.train = Random.modules\$Enrichment_Random_Modules[common_modules_random, ]
Random.modules.test = Random.modules.test[common_modules_random, ]
Random.model = calibrate_model(Input = Random.modules.train,
                               modules_metadata = train.meta,
                               label = "label", sampleIDs = "sampleIDs",
                               Metric = "Exact_match",
                               n_folds = folds, repetitions = 1)
Random.model[["ConfusionMatrix_training"]]

# validate on test set
Input.test <- t(Random.modules.test) %>% scale(center = TRUE, scale = TRUE)
Input.test<-Input.test[ , apply(Input.test, 2, function(x) !any(is.na(x)))]
Prediction.model <- predict(Random.model\$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing_random <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test.meta\$label))
Random.model\$ConfusionMatrix_testing_random <- ConfusionMatrix_testing_random
ConfusionMatrix_testing_random

#-------------------------PLOT PERFOMANCES-------------------------
train_name <- paste0(split_training * 100, "%") # MALTE TODO: maybe also escape
test_name <- paste0((1-split_training) * 100, "%")
metrics_plot <- plot_performance(trained_model = trained.model,
                                random_model = Random.model,
                                central_genes_model = NULL,
                                training_dataset_name = train_name,
                                testing_dataset_name = test_name,
                                subtypes = as.factor(train.meta\$label))
# save metrics plot
ggsave(file.path(PLOT_DIR, "metrics.png"), metrics_plot,
       width = 7.25, height = 5.25, dpi = 1200)

#-------------------------PLOT CLASS DISTRIBUTION------------------
density_plot_train <- plot_density_scores(trained_model = trained.model,
                                          spongEffects = train.modules,
                                          meta_data = train.meta,
                                          label = "label",
                                          sampleIDs = "sampleIDs")
# save class plot
ggsave(file.path(PLOT_DIR, "classification.png"), density_plot_train,
       width = 12, height = 8, dpi = 300)

#-------------------------PLOT TOP RESULTS-------------------------
lollipop_plot <- plot_modules(trained_model = trained.model, k_modules_red = 2,
                             k_modules = 7, text_size = 20)
# save lollipop plot
ggsave(file.path(PLOT_DIR, "lollipop.png"), lollipop_plot,
       width = 7.25, height = 5.25, dpi = 1200)

#-------------------------PLOT NETWORK FOR CENTRALITIES------------
n = 25
central_players <- central_gene_modules %>%
  arrange(desc(Weighted_Degree)) %>%
  dplyr::select(gene) %>% slice_head(n = n) %>% unlist
# filter network for central players, filter for circRNA-mRNA connections
network <- circ.mRNA.subnetwork(filtered_network_centralities\$Sponge.filtered, "circ")
subnetwork <- network %>%
  dplyr::filter(geneA %in% central_players | geneB %in% central_players)
network_plot <- plot_network(subnetwork,
                             annotation = gene.ens.all,
                             node_label_size = 40, edge_label_size = 0)
# save network
visNetwork::visSave(network_plot, file = file.path(PLOT_DIR, paste0("top_", n, "_ceRNA_centralities.html")))

#-------------------------EXTRACT TOP circRNAs---------------------
Variable.importance <- importance(trained.model\$Model\$finalModel) %>% as.data.frame() %>% 
  tibble::rownames_to_column("Module") %>% arrange(desc(MeanDecreaseGini))
k = 6
candidates <- Variable.importance[1:k, "Module"]
candidates <- gsub("`", "", candidates)
# candidate modules
candidate.modules <- Sponge.modules[candidates]
# continue with top 2 (highest gini indices)
sponged.genes <- candidate.modules[candidates]

#-------------------------PLOT TOP circRNA modules-----------------
top_network_plot <- plot_network(ceRNA_network = network %>%
                                dplyr::filter(geneA %in% candidates | geneB %in% candidates),
                                annotation = gene.ens.all,
                                node_label_size = 40, edge_label_size = 0)
# plot top k centralities network
visNetwork::visSave(top_network_plot,
                    file = file.path(PLOT_DIR, paste0("top_", k, "_ceRNA_centralities.html")))

#-------------------------MODULE EXPRESSION------------------------
# generate module expression for each candidate
candidate_plots <- list()
for (candidate in candidates) {
  module <- candidate.modules[candidate][[1]]
  module <- module[!grepl("c", module)]
  module_gene_plot <- plot_target_gene_expressions(candidate, module,
                                                  gene_expression = gene_expr,
                                                  meta = meta,
                                                  annotation = gene.ens.all,
                                                  unit = "counts",
                                                  log_transform = TRUE,
                                                  show_rows = FALSE)
  module_miRNA_plot <- plot_miRNAs_per_gene(candidate,
                                            genes_miRNA_interactions_manual,
                                            mir_expr = mir_expr, meta = meta,
                                            log_transform = TRUE)
  candidate_plots[[candidate]] <- ggarrange(module_gene_plot, module_miRNA_plot,
                                            ncol = 2, nrow = 1)
}
# arrange candidate plots
candidate_plots_arranged <- ggarrange(plotlist = candidate_plots,
                                      ncol = 1, nrow = k)
# save candidate heat maps
ggsave(file.path(PLOT_DIR, "candidates_hms.png"), candidate_plots_arranged,
       width = 5.25, height = 10.25, dpi = 1200)
       
#-------------------------spongEffects HMAP------------------------
train.hmap <- plot_hmap(trained_model = trained.model,
                        spongEffects = train.modules,
                        meta_data = train.meta,
                        label = "label", sampleIDs = "sampleIDs",
                        Modules_to_Plot = 2,
                        show.rownames = TRUE,
                        show.colnames = FALSE)
test.hmap <- plot_hmap(trained_model = trained.model,
                       spongEffects = test.modules,
                       meta_data = test.meta,
                       label = "label",
                       sampleIDs = "sampleIDs",
                       Modules_to_Plot = 2,
                       show.rownames = TRUE,
                       show.colnames = FALSE)
# save train and test heat maps
ggsave(file.path(PLOT_DIR, "train_hm.png"), train.hmap,
       width = 7.25, height = 5.25, dpi = 1200)
ggsave(file.path(PLOT_DIR, "test_hm.png"), test.hmap,
       width = 7.25, height = 5.25, dpi = 1200)

#-------------------------SAVE R-SESSION----------------------------
save.image(file.path(ROOT, "spongEffects.RData"))


################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

doparallel.version <- as.character(packageVersion('doParallel'))
dplyr.version <- as.character(packageVersion('dplyr'))
foreach.version <- as.character(packageVersion('foreach'))
ggplot2.version <- as.character(packageVersion('ggplot2'))
ggpubr.version <- as.character(packageVersion('ggpubr'))
gsva.version <- as.character(packageVersion('GSVA'))
pheatmap.version <- as.character(packageVersion('pheatmap'))
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
reshape2.version <- as.character(packageVersion('reshape2'))
rtracklayer.version <- as.character(packageVersion('rtracklayer'))
sponge.version <- as.character(packageVersion('SPONGE'))
visnetwork.version <- as.character(packageVersion('visNetwork'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-gsva:', gsva.version)
        paste('    bioconductor-rtracklayer', rtracklayer.version)
        paste('    bioconductor-sponge', sponge.version)
        paste('    r-base:', r.version)
        paste('    r-doparallel', doparallel.version)
        paste('    r-dplyr', dplyr.version)
        paste('    r-foreach', foreach.version)
        paste('    r-ggplot2', ggplot2.version)
        paste('    r-ggpubr', ggpubr.version)
        paste('    r-pheatmap', pheatmap.version)
        paste('    r-reshape2', reshape2.version)
        paste('    r-visnetwork', visnetwork.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
