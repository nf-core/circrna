process FIND_CIRC_ANCHORS {
    tag "$meta.id"
    label "process_high"

    conda "bioconda::find_circ=1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/find_circ%3A1.2--hdfd78af_0' :
        'biocontainers/find_circ:1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}_anchors.qfa.gz"), emit: anchors
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2'
    """
    unmapped2anchors.py $bam | gzip > ${prefix}_anchors.qfa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_circ: $VERSION
    END_VERSIONS
    """
}
