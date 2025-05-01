process PSIRC_FLI {
    tag "${meta.id}"
    label 'process_high'

    container 'docker.io/nicotru/psirc'

    input:
    tuple val(meta), path(reads), path(bsj, stageAs: 'bsj_output')
    tuple val(meta2), path(transcriptome)

    output:
    //tuple val(meta), path("output/candidate_circ_junctions.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    VERSION = "1.0.0"
    """
    cp -rL ${bsj} output
    psirc -s -t ${task.cpus} ${transcriptome} output ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        psirc: ${VERSION}
    END_VERSIONS
    """
}
