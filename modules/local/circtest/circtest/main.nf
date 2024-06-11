process CIRCTEST_CIRCTEST {
    label 'process_medium'

    conda "conda-forge::r-base=4.2.2 conda-forge::r-aod=1.3.2 conda-forge::r-ggplot2=3.4.0 conda-forge::r-plyr=1.8.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' :
        'biocontainers/mulled-v2-c79b00aa4647c739dbe7e8480789d3ba67988f2e:0' }"

    input:
    tuple val(meta) , path(circ_counts)
    tuple val(meta2), path(gene_counts)
    tuple val(meta3), path(phenotype)

    output:
    path "*"           , emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'circtest.R'
}
