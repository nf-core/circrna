process DESEQ2_NORMALIZATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0cde4f448f7fb2a5421bb9bf747e2c17fc7413916790dacca37930691d8e01d/data'
        : 'community.wave.seqera.io/library/deseq2_normalization:d76beacafc69a865'}"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.id}.normalized_counts.tsv"), emit: normalized
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('deseq_normalization.R')

    stub:
    """
    touch ${meta.id}.normalized_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "\$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")"
        bioconductor-deseq2: "\$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")"
    END_VERSIONS

    """
}
