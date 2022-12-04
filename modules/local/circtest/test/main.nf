process CIRCTEST {
    label 'process_medium'

    conda (params.enable_conda ? "r-base r-aod r-ggplot2 r-plyr" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' :
        'quay.io/biocontainers/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' }"

    input:
    path(circ_csv)
    path(linear_csv)
    path(phenotype)

    output:
    path "*"           , emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    circ_test.R $circ_csv $linear_csv $phenotype

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        aod: \$(Rscript -e "library(aod); cat(as.character(packageVersion('aod')))")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        plyr: \$(Rscript -e "library(plyr); cat(as.character(packageVersion('plyr')))")
    END_VERSIONS
    """
}
