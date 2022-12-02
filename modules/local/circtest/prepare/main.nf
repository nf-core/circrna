process PREPARE_CLR_TEST {
    label 'process_medium'

    conda (params.enable_conda ? "r-base r-aod r-ggplot2 r-plyr" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' :
        'quay.io/biocontainers/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' }"

    input:
    path(gene_matrix)
    path(circrna_matrix)
    path(circ_host_map)
    path(gtf)

    output:
    path "circ.csv"    , emit: circular
    path "linear.csv"  , emit: linear
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_circ_test.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
