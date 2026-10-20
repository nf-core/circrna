process PSIRC_FLI {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5afa9a01488ef93da1c7d72dfb0f0c9d885e5933aa67aa4ccb52f3c98eeb0c/data'
        : 'community.wave.seqera.io/library/psirc:1.0.0--4c24a8f727283097'}"

    input:
    tuple val(meta), path(reads), path(bsj, stageAs: 'bsj_output')
    tuple val(meta2), path(transcriptome)

    output:
    tuple val(meta), path("${prefix}/full_length_isoforms.fa"), emit: fasta
    tuple val(meta), path("${prefix}/full_length_isoforms_alt_fsj_supporting_reads.sam"), emit: sam
    tuple val(meta), path("${prefix}/full_length_isoforms.tsv"), emit: tsv
    tuple val("${task.process}"), val('psirc'), val('1.0.0'), emit: versions_psirc, topic: versions
    tuple val("${task.process}"), val('kallisto_b'), eval('kallisto_b version | sed "s/.*version //"'), emit: versions_kallisto_b, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -rL ${bsj} ${prefix}
    psirc -s -t ${task.cpus} ${prefix} ${transcriptome} "\$(command -v kallisto_b)" ${reads}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/full_length_isoforms.fa
    touch ${prefix}/full_length_isoforms_alt_fsj_supporting_reads.sam
    touch ${prefix}/full_length_isoforms.tsv
    """
}
