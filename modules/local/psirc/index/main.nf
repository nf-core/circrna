process PSIRC_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    container "registry.hub.docker.com/bigdatainbiomedicine/psirc"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("psirc.index")

    script:
    """
    psirc-quant index -i psirc.index --make-unique $fasta
    """
}