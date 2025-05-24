process CIRCTOOLS_DETECT {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/circtools:2.0--pyhdfd78af_0'
        : 'biocontainers/circtools:2.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(paired, stageAs: 'paired.junctions'), path(mate1, stageAs: 'mate1.junctions'), path(mate2, stageAs: 'mate2.junctions')
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}_reads.junctions"), emit: reads
    tuple val(meta), path("${prefix}_coordinates.tsv"), emit: coordinates
    tuple val(meta), path("${prefix}_counts.tsv")     , emit: counts

    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = meta.strandedness ?: 'auto'
    def strand_args = strandedness == 'auto' || strandedness == 'unstranded' ? '-N' : strandedness == 'forward' ? '' : '-ss'

    def matefile_commands = meta.single_end
        ? ''
        : "printf ${mate1} > mate1file && printf ${mate2} > mate2file"

    def mate_args = meta.single_end ? '' : '-mt1 @mate1file -mt2 @mate2file -Pi'
    """
    printf "${paired}" > samplesheet
    ${matefile_commands}

    circtools detect @samplesheet ${mate_args} -D -an ${gtf} ${args} -F -M -k -Nr 1 1 -A ${fasta} ${strand_args} -T ${task.cpus}

    mv _tmp_circtools/tmp_printcirclines.[0-9A-Z][0-9A-Z][0-9A-Z][0-9A-Z][0-9A-Z][0-9A-Z] ${prefix}_reads.junctions
    mv CircCoordinates ${prefix}_coordinates.tsv
    mv CircRNACount ${prefix}_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circtools: \$(circtools -V)
    END_VERSIONS
    """
}
