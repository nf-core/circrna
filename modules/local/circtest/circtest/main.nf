process CIRCTEST_CIRCTEST {
    label 'process_medium'
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1a6c264440ced99bec3c70893f1738b90421307432663d64815ea4384a0a9856/data'
        : 'community.wave.seqera.io/library/circtest_circtest:ff686a5b50cb6d65'}"

    input:
    tuple val(meta) , path(gene_counts), path(circ_counts)
    tuple val(meta2), path(phenotype)

    output:
    tuple val(meta), path("${prefix}_summary.txt"), emit: summary
    tuple val(meta), path("*.pdf")                , emit: plots
    tuple val("${task.process}"), val('r-base'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('aod'), eval('Rscript -e "cat(as.character(packageVersion(\'aod\')))"'), emit: versions_aod, topic: versions
    tuple val("${task.process}"), val('plyr'), eval('Rscript -e "cat(as.character(packageVersion(\'plyr\')))"'), emit: versions_plyr, topic: versions
    tuple val("${task.process}"), val('ggplot2'), eval('Rscript -e "cat(as.character(packageVersion(\'ggplot2\')))"'), emit: versions_ggplot2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'circtest.sh'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_summary.txt
    touch circ_linear_ratio_plots.pdf
    touch circ_linear_line_plots.pdf

    """
}
