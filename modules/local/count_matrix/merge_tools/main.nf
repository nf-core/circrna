process MERGE_TOOLS {
    tag "$meta.id"
    label 'process_low'

    conda "r-base=3.6.3 python=2.7.15 r-argparser=0.6 r-dplyr=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5fbffedf7f529cf3c5093b976deb4290f5e1267a:3456f1432b1c9dad42815275abe2d6cb6f26fd94-0' :
        'quay.io/biocontainers/mulled-v2-5fbffedf7f529cf3c5093b976deb4290f5e1267a:3456f1432b1c9dad42815275abe2d6cb6f26fd94-0' }"

    input:
    tuple val(meta), path(bed)
    val(tool_filter)
    val(duplicates_fun)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## make list of files for R to read
    ls *.bed > samples.csv

    ## Add catch for empty bed file and delete
    bash ${workflow.projectDir}/bin/check_empty.sh

    ## Use intersection of "n" (params.tool_filter) circRNAs called by tools
    ## remove duplicate IDs, keep highest count.
    Rscript ${workflow.projectDir}/bin/consolidate_algorithms_intersection.R samples.csv $tool_filter $duplicates_fun
    mv combined_counts.bed ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        argparser: \$(Rscript -e "library(arparser); cat(as.character(packageVersion('argparser')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
