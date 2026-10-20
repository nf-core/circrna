process CIRCTEST_CIRCTEST {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1a6c264440ced99bec3c70893f1738b90421307432663d64815ea4384a0a9856/data'
        : 'community.wave.seqera.io/library/circtest_circtest:ff686a5b50cb6d65'}"

    input:
    tuple val(meta), path(gene_counts), path(circ_counts)
    tuple val(meta2), path(phenotype)

    output:
    tuple val(meta), path("${prefix}_summary.txt"), emit: summary
    tuple val(meta), path("*.pdf"), emit: plots
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template('circtest.R')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_summary.txt
    touch circ_linear_ratio_plots.pdf
    touch circ_linear_line_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "\$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")"
        aod: "\$(Rscript -e "cat(as.character(packageVersion('aod')))")"
        plyr: "\$(Rscript -e "cat(as.character(packageVersion('plyr')))")"
        ggplot2: "\$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")"
    END_VERSIONS

    """
}
