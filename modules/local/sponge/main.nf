process SPONGE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-sponge_r-doparallel_r-foreach_r-visnetwork:379bbce95f2c4913' :
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-sponge_r-doparallel_r-foreach_r-visnetwork:f38bbe7441604bcc' }"

    input:
    tuple val(meta), path(binding_sites)
    tuple val(meta), path(gene_expr)
    tuple val(meta), path(mirna_expr)


    output:
    path "versions.yml"                             , emit: versions
    tuple val(meta), path("sponge.RData")           , emit: sponge_data
    path "circRNAs_as_ceRNAs.tsv"                   , emit: cernas
    path "total_plots.pdf"                          , emit: total_plots
    path "simulation_plots.pdf"                     , emit: simulation_plots
    path "circRNA_plots.pdf"                        , emit: circ_plots


    when:
    task.ext.when == null || task.ext.when

    script:
    template 'sponge.R'

    stub:
    """
    touch "sponge.RData"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        bioconductor-sponge: \$(Rscript -e "library(SPONGE); cat(as.character(packageVersion('SPONGE')))")
        r-visnetwork: \$(Rscript -e "library(visNetwork); cat(as.character(packageVersion('visNetwork')))")
        r-doparallel: \$(Rscript -e "library(doParallel); cat(as.character(packageVersion('doParallel')))")
        r-foreach: \$(Rscript -e "library(foreach); cat(as.character(packageVersion('foreach')))")
    END_VERSIONS
    """
}
