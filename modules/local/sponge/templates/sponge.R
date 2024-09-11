#!/usr/bin/env Rscript

library(SPONGE)
library(visNetwork)
library(doParallel)
library(foreach)

circ_mRNA_subnetwork <- function(interactions, pattern) {
    return(interactions[
        (grepl(pattern, interactions\$geneA) & !grepl(pattern, interactions\$geneB)) |
        (grepl(pattern, interactions\$geneB) & !grepl(pattern, interactions\$geneA)),
    ])
}

# load binding sites matrix
bindingsites <- read.table("${binding_sites}", sep = "\\t", header = TRUE, check.names = FALSE)
rownames(bindingsites) <- bindingsites[, 1]
bindingsites <- bindingsites[, -1]
bindingsites_matrix <- as.matrix(bindingsites)


# load gene expression
genes <- read.csv("${gene_expr}", sep = "\\t", row.names = "tx")
genes <- subset(genes, select = -gene_id)
gene_expr <- t(as.matrix(genes))


# load mirna expression
mirna <- read.csv("${mirna_expr}", sep = "\\t", row.names = "miRNA")
mirna[is.na(mirna)] <- 0
mir_expr <- t(as.matrix(mirna))


# ensure that both matrices have the same samples
common_samples <- intersect(rownames(gene_expr), rownames(mir_expr))

# subset and reorder the rows (samples) of both matrices to match each other (as specified in SPONGE docs)
gene_expr <- gene_expr[common_samples, ]
mir_expr <- mir_expr[common_samples, ]


######################################
#                A                   #
######################################
logging.file <- "sponge.log"
# number of cores from args
cores <- ${params.sponge_cpus}

cl <- makeCluster(cores, outfile=logging.file)
registerDoParallel(cl)

gene_expriRNA_candidates <- SPONGE::sponge_gene_miRNA_interaction_filter(
    gene_expr = gene_expr,
    mir_expr = mir_expr,
    mir_predicted_targets = bindingsites_matrix,
    F.test = FALSE,
    coefficient.threshold = -0.01,
    coefficient.direction = NULL
)


    # F.test.p.adj.threshold = ${params.sponge_f_test_pval},
    # coefficient.threshold = ${params.sponge_coeff_threshold},
    # elastic.net = ("${params.sponge_elastic_net}" == "true")

######################################
#                B                   #
######################################
ceRNA_interactions <- SPONGE::sponge(
    gene_expr = gene_expr,
    mir_expr = mir_expr,
    mir_interactions = gene_expriRNA_candidates
)


######################################
#                C                   #
######################################

mscor_null_model <- SPONGE::sponge_build_null_model(
    number_of_datasets = 100,
    number_of_samples = nrow(gene_expr)
)

sim_plot <- SPONGE::sponge_plot_simulation_results(mscor_null_model)
pdf("simulation_plots.pdf")
plot(sim_plot)
dev.off()

ceRNA_interactions_sign <- SPONGE::sponge_compute_p_values(
    sponge_result = ceRNA_interactions,
    null_model = mscor_null_model
)

######################################
#                D                   #
######################################
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign\$p.adj < ${params.sponge_fdr}), ]
min.interactions <- 5000

if (nrow(ceRNA_interactions_fdr) < min.interactions && nrow(ceRNA_interactions_sign) > min.interactions) {
    print("Warning: fdr setting too strict, no significant interactions detected; min of padj is:")
    print(min(ceRNA_interactions_sign\$p.adj))
    print("adjusting...")
    fdr <- min(ceRNA_interactions_sign\$p.adj) * 1.01
    ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign\$p.adj < fdr), ]
    while(nrow(ceRNA_interactions_fdr) < min.interactions) {
        fdr <- fdr * 1.01
        ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign\$p.adj < fdr), ]
    }
    cat("adjusted fdr to :", fdr, "to allow for a minimum", min.interactions, "interactions", "\\n")
    cat("current fdr relevant interactions:", nrow(ceRNA_interactions_fdr))
    ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr\$p.adj), ]
}

# MOST SIGNIFICANT SPONGES
network_centralities <- SPONGE::sponge_node_centralities(ceRNA_interactions_fdr)
ceRNA_interactions_weight <- ceRNA_interactions_fdr
ceRNA_interactions_weight\$weight <- -log10(ceRNA_interactions_fdr\$p.adj)
weighted_network_centralities <- SPONGE::sponge_node_centralities(ceRNA_interactions_weight)
weighed_network_plot <- SPONGE::sponge_plot_network_centralities(weighted_network_centralities, top = 3)
pdf("total_plots.pdf")
plot(weighed_network_plot)
dev.off()

# CIRC RNA MRNA ONLY
circ_mRNA_only <- circ_mRNA_subnetwork(ceRNA_interactions_fdr, "circ")
write.table(circ_mRNA_only, "circRNAs_as_ceRNAs.tsv", sep = "\\t", row.names = FALSE)

# TODO: REMOVE
write.table(ceRNA_interactions, "ceRNAs_interactions.tsv", sep = "\\t", row.names = FALSE)

# NETWORK ANALYSIS CIRC
ceRNA_interactions_circ_weight <- circ_mRNA_only
ceRNA_interactions_circ_weight\$weight <- -log10(circ_mRNA_only\$p.val)
weighted_network_centralities_circ <- SPONGE::sponge_node_centralities(ceRNA_interactions_circ_weight)

# plot top n samples
n = 10 # TODO: param

# betweeness
top_network_plot_btw <- SPONGE::sponge_plot_network_centralities(
    weighted_network_centralities_circ,
    top = n,
    measure = "btw"
)

# eigenvector
top_network_plot_ev <- SPONGE::sponge_plot_network_centralities(
    weighted_network_centralities_circ,
    top = n,
    measure = "ev"
)

# counts
top_network_plot_c <- SPONGE::sponge_plot_network_centralities(
    weighted_network_centralities_circ,
    top = n,
    measure = "count"
)

pdf("circRNA_plots.pdf")
plot(top_network_plot_btw)
plot(top_network_plot_ev)
plot(top_network_plot_c)
dev.off()
stopCluster(cl)

# save R objects for sponge effects
save.image(file =  "sponge.RData")

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

sponge.version <- as.character(packageVersion('SPONGE'))
visnetwork.version <- as.character(packageVersion('visNetwork'))
doparallel.version <- as.character(packageVersion('doParallel'))
foreach.version <- as.character(packageVersion('foreach'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-sponge:', sponge.version),
        paste('    r-visnetwork:', visnetwork.version),
        paste('    r-doparallel:', doparallel.version),
        paste('    r-foreach:', foreach.version)
    ),
'versions.yml')


################################################
################################################
################################################
################################################
################################################
