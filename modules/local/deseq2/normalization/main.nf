process DESEQ2_NORMALIZATION {
    tag "$meta.id"
    label 'process_single'
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0cde4f448f7fb2a5421bb9bf747e2c17fc7413916790dacca37930691d8e01d/data'
        : 'community.wave.seqera.io/library/deseq2_normalization:d76beacafc69a865'}"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.id}.normalized_counts.tsv"), emit: normalized
    tuple val("${task.process}"), val('r-base'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval('Rscript -e "cat(as.character(packageVersion(\'DESeq2\')))"'), emit: versions_bioconductor_deseq2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq_normalization.sh'

    stub:
    """
    touch ${meta.id}.normalized_counts.tsv

    """
}
