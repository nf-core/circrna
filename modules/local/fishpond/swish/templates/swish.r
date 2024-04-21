#!/usr/bin/env Rscript --vanilla

library(fishpond)

se <- readRDS('$experiment')

se <- scaleInfReps(se)
se <- labelKeep(se)
se <- swish(se, x="$column")

saveRDS(se, '${meta.id}.rds')

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-fishpond:', packageVersion('fishpond'))
    ),
'versions.yml')
