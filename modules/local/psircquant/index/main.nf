process PSIRCQUANT_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5afa9a01488ef93da1c7d72dfb0f0c9d885e5933aa67aa4ccb52f3c98eeb0c/data'
        : 'community.wave.seqera.io/library/psirc:1.0.0--4c24a8f727283097'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("psirc.index"), emit: index
    tuple val("${task.process}"), val('psirc-quant'), eval('psirc-quant version | sed -n "s/^psirc-quant, version \\([0-9.]*\\).*$/\\1/p"'), emit: versions_psircquant, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    psirc-quant index -i psirc.index --make-unique ${fasta}

    """

    stub:
    """
    touch psirc.index

    """
}
