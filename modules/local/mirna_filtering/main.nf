process MIRNA_FILTERING {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed8bfdd256acb1cb60f9b1252aafd0dd726fd61da2b62669e78741b744ed7403/data'
        : 'community.wave.seqera.io/library/mirna_filtering:6d818378228715c0'}"

    input:
    tuple val(meta), path(normalized_counts)
    val mirna_min_sample_percentage
    val mirna_min_reads

    output:
    tuple val(meta), path("${meta.id}.normalized_counts_filtered.tsv"), emit: filtered
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('mirna_filtering.R')

    stub:
    """
    touch ${meta.id}.normalized_counts_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "\$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")"
    END_VERSIONS

    """
}
