process PSIRC_TRANSCRIPTOME {
    tag "${meta.id}"
    label 'process_single'

    container "registry.hub.docker.com/bigdatainbiomedicine/psirc"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.transcriptome.fasta")

    script:
    """
    psirc_transcriptome.pl ${fasta} ${gtf} ${meta.id}.transcriptome.fasta
    """
}