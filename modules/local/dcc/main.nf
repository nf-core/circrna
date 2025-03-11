process DCC {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::circtools=2.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/circtools:2.0--pyhdfd78af_0'
        : 'biocontainers/circtools:2.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(paired), path(mate1), path(mate2)
    path fasta
    path gtf

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = meta.strandedness ?: 'auto'
    def strand_args = strandedness == 'auto' || strandedness == 'unstranded' ? '-N' : strandedness == 'forward' ? '' : '-ss'

    def matefile_commands = meta.single_end ? '' :
        "printf ${mate1} > mate1file && printf ${mate2} > mate2file"

    def mate_args = meta.single_end ? '' : '-mt1 @mate1file -mt2 @mate2file -Pi'

    """
    printf "${paired}" > samplesheet
    ${matefile_commands}

    circtools detect @samplesheet ${mate_args} -D -an ${gtf} ${args} -F -M -Nr 1 1 -A ${fasta} ${strand_args} -T ${task.cpus}

    awk '{print \$6}' CircCoordinates >> strand
    paste CircRNACount strand | tail -n +2 | awk -v OFS="\\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circtools: \$(circtools -V)
    END_VERSIONS
    """
}
