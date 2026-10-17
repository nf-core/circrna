process CIRCTEST_CIRCTEST {
    label 'process_medium'

    conda "conda-forge::r-base=4.2.2 conda-forge::r-aod=1.3.2 conda-forge::r-ggplot2=3.4.0 conda-forge::r-plyr=1.8.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' :
        'biocontainers/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' }"

    input:
    tuple val(meta) , path(gene_counts), path(circ_counts)
    tuple val(meta2), path(phenotype)

    output:
    tuple val(meta), path("${prefix}_summary.txt"), emit: summary
    tuple val(meta), path("*.pdf")                , emit: plots
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'circtest.R'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_summary.txt
    touch circ_linear_ratio_plots.pdf
    touch circ_linear_line_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | cut -d' ' -f3)
        aod: \$(Rscript -e 'cat(as.character(packageVersion("aod")))')
        plyr: \$(Rscript -e 'cat(as.character(packageVersion("plyr")))')
        ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))')
    END_VERSIONS
    """
}
