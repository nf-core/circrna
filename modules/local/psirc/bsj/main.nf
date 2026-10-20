process PSIRC_BSJ {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5afa9a01488ef93da1c7d72dfb0f0c9d885e5933aa67aa4ccb52f3c98eeb0c/data'
        : 'community.wave.seqera.io/library/psirc:1.0.0--4c24a8f727283097'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(transcriptome), path(index)

    output:
    tuple val(meta), path("output"), emit: output
    tuple val(meta), path("output/candidate_circ_junctions.bed"), emit: bed
    tuple val("${task.process}"), val('psirc'), val('1.0.0'), emit: versions_psirc, topic: versions
    tuple val("${task.process}"), val('kallisto_b'), eval('kallisto_b version | sed "s/.*version //"'), emit: versions_kallisto_b, topic: versions

    script:
    """
    mkdir -p output
    psirc -o output -t ${task.cpus} ${transcriptome} "\$(command -v kallisto_b)" ${reads}
    """

    stub:
    """
    mkdir -p output
    touch output/candidate_circ_junctions.bed
    """
}
