process PSIRC_INDEX {
    tag "${meta.id}"
    label 'process_high'

    container "registry.hub.docker.com/bigdatainbiomedicine/psirc"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("psirc.index"), emit: index
    path "versions.yml",                  emit: versions

    script:
    def VERSION = '1.0'
    """
    psirc-quant index -i psirc.index --make-unique $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc-quant: $VERSION
    END_VERSIONS
    """
}
