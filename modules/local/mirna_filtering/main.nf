process MIRNA_FILTERING {
    tag "${meta.id}"
    label 'process_single'
    shell '/bin/bash'

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
    tuple val("${task.process}"), val('r-base'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('mirna_filtering.sh')

    stub:
    """
    touch ${meta.id}.normalized_counts_filtered.tsv

    """
}
