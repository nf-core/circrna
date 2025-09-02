process PSIRC_BSJ {
    tag "${meta.id}"
    label 'process_high'

    container 'docker.io/nicotru/psirc'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(transcriptome), path(index)

    output:
    tuple val(meta), path("output"), emit: output
    tuple val(meta), path("output/candidate_circ_junctions.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    VERSION = "1.0.0"
    """
    mkdir -p output
    psirc -o output -t ${task.cpus} ${transcriptome} ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc: ${VERSION}
    END_VERSIONS
    """
}
