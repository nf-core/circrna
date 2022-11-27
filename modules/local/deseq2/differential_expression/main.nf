process DESEQ2_DIFFERENTIAL_EXPRESSION {
    label 'process_medium'

    conda (params.enable_conda ? "r-base=3.6.3 conda-forge::r-argparser=0.6 conda-forge::r-dplyr=1.0.5 conda-forge::r-ggplot2=3.3.3 r-ggpubr=0.4.0 conda-forge::r-gplots=3.1.1 conda-forge::r-pheatmap=1.0.12 r-plyr=1.8.6 r-pvclust=2.2_0 r-rcolorbrewer=1.1_2 conda-forge::r-circlize=0.4.12 bioconductor-biomart=2.42.0 bioconductor-complexheatmap=2.2.0 bioconductor-deseq2=1.26.0 bioconductor-enhancedvolcano=1.4.0 bioconductor-ihw=1.14.0 bioconductor-org.hs.eg.db=3.10.0 bioconductor-pcatools=1.2.0 bioconductor-tximport=1.14.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/biocontainers/mulled-v2-04b2ef814c9c6ab8c196c3e372521b88160dc260:e0cb4046baee3fd35fdbf883ba8af34e3e8af2e8-0' :
        'quay.io/biocontainers/mulled-v2-04b2ef814c9c6ab8c196c3e372521b88160dc260:e0cb4046baee3fd35fdbf883ba8af34e3e8af2e8-0' }"

    input:
    file(gene_matrix)
    file(phenotype)
    file(circrna_matrix)
    val(species)

    output:
    path "circRNA" , emit: circular_results
    path "RNA-Seq"       , emit: linear_results
    path "boxplots"             , emit: boxplots
    path "DESeq2_QC"    , emit: qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## prepDE && circRNA counts headers are sorted where uppercase preceedes lowercase i.e Z before a
    ## reformat the phenotype file to match the order of the samples.
    head -n 1 $phenotype > header
    tail -n +2 $phenotype | LC_COLLATE=C sort > sorted_pheno
    cat header sorted_pheno > tmp && rm phenotype.csv && mv tmp phenotype.csv

    Rscript DEA.R $gene_matrix $phenotype $circrna_matrix $species ensembl_database_map.txt
    """
}
