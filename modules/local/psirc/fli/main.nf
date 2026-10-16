process PSIRC_FLI {
    tag "${meta.id}"
    label 'process_high'

    container 'docker.io/nicotru/psirc'

    input:
    tuple val(meta), path(reads), path(bsj, stageAs: 'bsj_output')
    tuple val(meta2), path(transcriptome)

    output:
    tuple val(meta), path("${prefix}/full_length_isoforms.fa"), emit: fasta
    tuple val(meta), path("${prefix}/full_length_isoforms_alt_fsj_supporting_reads.sam"), emit: sam
    tuple val(meta), path("${prefix}/full_length_isoforms.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    VERSION = "1.0.0"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -rL ${bsj} ${prefix}
    psirc -s -t ${task.cpus} ${transcriptome} ${prefix} ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc: ${VERSION}
    END_VERSIONS
    """
}
