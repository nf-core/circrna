process SPONGE_EFFECTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-gsva_bioconductor-rtracklayer_bioconductor-sponge_r-base_pruned:ac5be0e145e964fc' :
        'community.wave.seqera.io/library/bioconductor-gsva_bioconductor-rtracklayer_bioconductor-sponge_r-base_pruned:5fad30741cf73551' }"



    input:
    tuple val(meta), path(sponge_data)

    output:
    // TODO


    when:
    task.ext.when == null || task.ext.when

    script:
    template 'sponge_effects.R'

    stub:
    """
    touch "test.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtracklayer: \$(Rscript -e "library(foreach); cat(as.character(packageVersion('foreach')))")
        bioconductor-sponge: \$(Rscript -e "library(SPONGE); cat(as.character(packageVersion('SPONGE')))")
        r-doparallel: \$(Rscript -e "library(doParallel); cat(as.character(packageVersion('doParallel')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-foreach: \$(Rscript -e "library(foreach); cat(as.character(packageVersion('foreach')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-ggpubr: \$(Rscript -e "library(ggpubr); cat(as.character(packageVersion('ggpubr')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
        r-reshape2: \$(Rscript -e "library(reshape2); cat(as.character(packageVersion('reshape2')))")
        r-visnetwork: \$(Rscript -e "library(visNetwork); cat(as.character(packageVersion('visNetwork')))")
    END_VERSIONS
    """
}
