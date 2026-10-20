process COMPUTE_CORRELATIONS {
    tag "$meta.id"
    label 'process_single'
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c3c7441de9df5430189904a4f577ac48cabbfbef61e64e90c1e401b28748ca1f/data'
        : 'community.wave.seqera.io/library/compute_correlations:b3dacbb39dbd06ae'}"

    input:
    tuple val(meta),  path(bindingsites)
    tuple val(meta2), path(mirna_expression)
    tuple val(meta3), path(transcript_rds)

    output:
    tuple val(meta), path("*.tsv"), emit: correlations, optional: true
    tuple val("${task.process}"), val('r-base'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('bioconductor-fishpond'), eval('Rscript -e "cat(as.character(packageVersion(\'fishpond\')))"'), emit: versions_bioconductor_fishpond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'compute_correlations.sh'

    stub:
    """
    touch ${meta.id}.circrna_correlation.tsv

    """
}
