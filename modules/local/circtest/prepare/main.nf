process CIRCTEST_PREPARE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/476ef7843285cc28e0778cb22a101abe14d3c218fb74bd881eb5fc22ceb750c9/data'
        : 'community.wave.seqera.io/library/circtest_prepare:35adc3cf2748fa7e'}"

    input:
    tuple val(meta), path(gene_counts), path(circ_counts)

    output:
    tuple val(meta), path('*_genes.tsv'), path('*_circs.tsv'), emit: counts, optional: true

    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template('prepare.R')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_genes.tsv
    touch ${prefix}_circs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "\$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")"
    END_VERSIONS

    """
}
