process SPONGE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-sponge_r-doparallel_r-foreach_r-visnetwork:086ecad7b1034ff0' :
        'community.wave.seqera.io/library/bioconductor-sponge_r-doparallel_r-foreach_r-visnetwork:e79cb4c59aecf7ba' }"

    input:
    tuple val(meta), path(binding_sites)
    tuple val(meta), path(gene_expr)
    tuple val(meta), path(mirna_expr)


    output:
    tuple val(meta), path("sponge.RData")           , emit: sponge_data
    path "circRNAs_as_ceRNAs.tsv"                   , emit: cernas
    path "total_plots.pdf"                          , emit: total_plots
    path "simulation_plots.pdf"                     , emit: simulation_plots
    path "circRNA_plots.pdf"                        , emit: circ_plots
    path "versions.yml"                             , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    template 'sponge.R'

    stub:
    """
    touch "sponge.RData"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-sponge: \$(Rscript -e "library(SPONGE); cat(as.character(packageVersion('SPONGE')))")
        r-visnetwork: \$(Rscript -e "library(visNetwork); cat(as.character(packageVersion('visNetwork')))")
        r-doparallel: \$(Rscript -e "library(doParallel); cat(as.character(packageVersion('doParallel')))")
        r-foreach: \$(Rscript -e "library(foreach); cat(as.character(packageVersion('foreach')))")
    END_VERSIONS
    """
}
