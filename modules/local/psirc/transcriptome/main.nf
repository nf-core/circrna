process PSIRC_TRANSCRIPTOME {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/nicotru/psirc'

    input:
    tuple val(meta) , path(gtf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("$prefix.$suffix"), emit: transcriptome
    path "versions.yml"                     , emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    suffix = task.ext.suffix ?: 'fasta'
    VERSION = "1.0.0"
    """
    create_custom_transcriptome_fa.pl $fasta $gtf $prefix.$suffix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc-quant: $VERSION
    END_VERSIONS
    """
}
