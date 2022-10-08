#!/usr/bin/env R

args = commandArgs(trailingOnly = TRUE)

circ = read.table(args[1], header=T, sep=",")
linear = read.table(args[2], header=T, sep=",")
pheno = read.table(args[3], header=T, sep=",", row.names = "Sample_ID")


# No need to enforce any filtering for circTest.
# 'filter.sample' - this has been applied to called circs using the 'tool.filter' param
# 'filter.count' - this has been applied to called circs using bsj_filter param
# 'percentage' - set to extremely low value (do not want to discard circRNAs - let the user inspect themselves).

# Need to apply the phenotype csv file correctly to circtest.
n_covars <- ncol(pheno)
if( n_covars == 2){
  covariate_1 <- as.factor(pheno[,1])
  covariate_2 <- as.factor(pheno[,2])
}else{
  covariate_1 <- as.factor(pheno[,1])
}

n_reps <- as.numeric(table(covariate_1)[1])

Circ_filtered <- Circ.filter(circ = circ, linear = linear, Nreplicates = n_reps, filter.sample = 1, filter.count = 1, percentage = 0.00001, circle_description = c(1:4))
Linear_filtered <- linear[rownames(Circ_filtered),]


# groups must be numerically encoded
group = as.numeric(covariate_1)
test <- Circ.test(Circ_filtered, Linear_filtered, group=group, circle_description = c(1:4))


# Apply pheno to output once more..

if( n_covars == 2 ){

  group_indicator1 <- as.character(covariate_1)
  group_indicator2 <- as.character(covariate_2)

  pdf("circ_linear_ratio_plots.pdf", width = 8, height = 10)
  for (i in rownames(test$summary_table))  {
    Circ.ratioplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=groupindicator1, groupindicator2 = group_indicator2,
  		  circle_description = c(1:4) )
  }
  dev.off()

  pdf("circ_linear_line_plots.pdf", width = 8, height = 10)
  for (i in rownames(test$summary_table))  {
  Circ.lineplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=group_indicator1, groupindicator2 = group_indicator2,
	circle_description = c(1:4) )
  }
  dev.off()

}else{

  group_indicator1 <- as.character(covariate_1)

  pdf("circ_linear_ratio_plots.pdf", width = 8, height = 10)
  for (i in rownames(test$summary_table))  {
    Circ.ratioplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=group_indicator1,
  		lab_legend = colnames(pheno)[1],  circle_description = c(1:4) )
  }
  dev.off()

  pdf("circ_linear_line_plots.pdf", width = 8, height = 10)
  for (i in rownames(test$summary_table))  {
  Circ.lineplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=group_indicator1,
	circle_description = c(1:4) )
  }
  dev.off()
}


# include variables, makes life easier in case user wishes to report bugs to workflow.
save.image("circ_test.RData")
