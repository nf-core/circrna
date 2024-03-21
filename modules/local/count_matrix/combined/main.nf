process COUNTS_COMBINED {
    label 'process_low'

    conda "conda-forge::r-base=3.6.3 conda-forge::python=2.7.15 conda-forge::r-argparser=0.6 conda-forge::r-dplyr=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5fbffedf7f529cf3c5093b976deb4290f5e1267a:3456f1432b1c9dad42815275abe2d6cb6f26fd94-0' :
        'quay.io/biocontainers/mulled-v2-5fbffedf7f529cf3c5093b976deb4290f5e1267a:3456f1432b1c9dad42815275abe2d6cb6f26fd94-0' }"

    input:
    path(bed)

    output:
    path("circRNA_matrix.txt"), emit: dea_matrix
    path("count_matrix.txt")  , emit: clr_matrix
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    circRNA_counts_matrix.py > matrix.txt
    ## handle non-canon chromosomes here (https://stackoverflow.com/questions/71479919/joining-columns-based-on-number-of-fields)
    n_samps=\$(ls *.bed | wc -l)
    canon=\$(awk -v a="\$n_samps" 'BEGIN {print a + 4}')
    awk -v n="\$canon" '{ for (i = 2; i <= NF - n + 1; ++i) { \$1 = \$1"-"\$i; \$i=""; } } 1' matrix.txt | awk -v OFS="\\t" '\$1=\$1' > circRNA_matrix.txt
    reformat_count_matrix.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | cut -d' ' -f3 | sed 's/,//g')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
