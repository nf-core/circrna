process PSIRC_INDEX {
    tag "${meta.id}"
    label 'process_high'

    container 'docker.io/nicotru/psirc'

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("*.all_possible_bsj_targets.fa.index"), path("*.FirstandLastExons_entities.fa.index"), emit: index
    path "versions.yml", emit: versions

    script:
    VERSION = "1.0.0"
    """
    psirc -i $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc: $VERSION
    END_VERSIONS
    """
}
