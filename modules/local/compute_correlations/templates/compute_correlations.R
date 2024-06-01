#!/usr/bin/env Rscript

dataset <- read.table('${params.input}', sep = ",", header=T, stringsAsFactors = F, check.names = F)
samples <- dataset\$sample

pairBindSites <- read.table('$bindingsites', header = F, sep = "\\t", stringsAsFactors = F, check.names = F)
pairBindSites <- as.data.frame.matrix(table(pairBindSites[,2], pairBindSites[,1]))

if (!nrow(pairBindSites) > 1) {  # check if there are samples that survived majority vote
  stop("\\nMajority vote resulted in no bindingsites! Check majority.tsv and adjust the parameter 'mirn_vote' in your nextflow.config to a lower number.")
}

miRNA_expression <- read.table('$filtered_mirnas', header = T, stringsAsFactors = F, check.names = F, sep = "\\t")


if(length(samples) < 5){
  stop("Cannot perform correlation on less than 5 samples")
}

circRNA_expression <- read.table('$filtered_circrnas', header = T, stringsAsFactors = F, check.names = F, sep="\\t")

circRNA_expression\$tx <- sub("^circ_", "", circRNA_expression\$tx)
rownames(circRNA_expression) <- circRNA_expression\$tx


# filter expression for miRNA binding pairs
valid_circRNAs <- intersect(circRNA_expression\$tx, rownames(pairBindSites))
valid_miRNAs <- intersect(miRNA_expression\$miRNA, colnames(pairBindSites))

circRNA_expression <- circRNA_expression[valid_circRNAs,]
miRNA_expression <- miRNA_expression[miRNA_expression\$miRNA %in% valid_miRNAs,]

# check for annotation
annotation <- "circBaseID" %in% colnames(circRNA_expression)  #TODO Do we need this?

header <- "circRNA\\tmiRNA\\tcircRNA_miRNA_ratio\\tmiRNA_binding_sites\\tpearson_R\\tcorr_pval\\tRSS_norm\\tintercept\\tintercept_pval\\tslope\\tslope_pval\\tadj_r_squared"
#write(header, file=paste0("filtered_circRNA_miRNA_correlation_libSizeEstNorm_directwritten.tsv"), append = F)

miRNA_for_row <- function(miRNA_expr_line, circRNA, circRNA_counts){
  
  mirna <- as.character(miRNA_expr_line[1])
  # get sample counts for current miRNA
  miRNA_counts <- miRNA_expr_line[-1]
    
  miRNA_counts <- data.frame(sample = as.character(names(miRNA_counts)), "miRNA_counts" = as.numeric(unname(miRNA_counts)))
  # compute circRNA expression vs. miRNA expression
  joined_counts <- merge(miRNA_counts, circRNA_counts, by="sample")

  # analyse circRNA/miRNA ratio
  mean_circRNA_counts <- mean(joined_counts\$circRNA_counts)
  mean_miRNA_counts <- mean(joined_counts\$miRNA_counts)
  circRNA_miRNA_ratio <- mean_circRNA_counts/mean_miRNA_counts
  
  binding_sites <- pairBindSites[circRNA, mirna]

  # compute circRNA-miRNA correlation for all samples
  cor_res <- cor.test(joined_counts\$miRNA_counts, joined_counts\$circRNA_counts,  method = "pearson", use = "complete.obs")
  corr_R <- as.numeric(as.character(cor_res\$estimate))
  corr_pval <- as.numeric(as.character(cor_res\$p.value))

  # compute linear regression
  regression_model <- lm(miRNA_counts~circRNA_counts, data = joined_counts)
  intercept <- summary(regression_model)\$coefficients[1,1]
  intercept_pval <- summary(regression_model)\$coefficients[1,4]
  slope <- summary(regression_model)\$coefficients[2,1]
  slope_pval <- summary(regression_model)\$coefficients[2,4]
  adj_r_squared <- summary(regression_model)\$adj.r.squared
  
  # compute residuals sum of squares
  # normalize counts for residuals sum of squares
  normalized_counts <- joined_counts[,c("circRNA_counts", "miRNA_counts")]
  min_circRNA_counts <- min(normalized_counts\$circRNA_counts)
  max_circRNA_counts <- max(normalized_counts\$circRNA_counts)
  normalized_counts[,"circRNA_counts"] <- (normalized_counts[,"circRNA_counts"] - min_circRNA_counts)/(max_circRNA_counts - min_circRNA_counts)
  min_miRNA_counts <- min(normalized_counts\$miRNA_counts)
  max_miRNA_counts <- max(normalized_counts\$miRNA_counts)
  normalized_counts[,"miRNA_counts"] <- (normalized_counts[,"miRNA_counts"] - min_miRNA_counts)/(max_miRNA_counts - min_miRNA_counts)
  norm_reg_model <- lm(miRNA_counts~circRNA_counts, data = normalized_counts)
  RSS_norm <- sum(norm_reg_model\$residuals^2)

  res <- data.frame(circRNA = as.character(circRNA), miRNA = as.character(mirna), 
                  circRNA_miRNA_ratio = as.numeric(circRNA_miRNA_ratio), 
                  miRNA_binding_sites = as.numeric(binding_sites), 
                  pearson_R = as.numeric(corr_R), corr_pval = as.numeric(corr_pval), 
                  RSS_norm = RSS_norm, intercept = intercept, 
                  intercept_pval = intercept_pval, slope = slope, 
                  slope_pval = slope_pval, adj_r_squared = adj_r_squared)
  
  return(res)
}

circRNA_for_row <- function(circRNA_expr_line){
  # get coordinations of current circRNA
  circRNA <- as.character(circRNA_expr_line[1])
  
  # annotate if possible
  if (annotation){
      an = circRNA_expression[circRNA,"circBaseID"]
      if (an != "None") {
            message("processing: ", circRNA, " (", an, ")")
            circRNA = an
      }
    
  } else {
      message("processing: ", circRNA)
  }
  
  # extract pure counts only
  circRNA_counts <- circRNA_expr_line[samples]  # passt
  
  # get sample counts for current circRNA
  circRNA_counts <- data.frame(sample = as.character(names(circRNA_counts)), "circRNA_counts" = as.numeric(unname(circRNA_counts)))
  
  if (all(circRNA_counts\$circRNA_counts != 0)) {
    res_list <- apply(miRNA_expression, 1, FUN = miRNA_for_row, circRNA, circRNA_counts)
    res_df <- do.call(rbind, res_list)
	return(res_df)
  }

}

correlations_list <- apply(circRNA_expression, MARGIN = 1, circRNA_for_row)
correlations_df <- do.call(rbind, correlations_list)
correlations_df\$adj_pval <- p.adjust(correlations_df\$corr_pval, method = "BH")

write.table(correlations_df, file=paste0("${meta.id}.circrna_correlation.tsv"), sep = "\\t", quote = F, row.names = F, col.names = T)


################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version)
    ), 
'versions.yml')
