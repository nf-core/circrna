process DCC {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::circtools=2.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/circtools:2.0--pyhdfd78af_0'
        : 'biocontainers/circtools:2.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(pairs), path(mate1), path(mate2)
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

    // Define variables based on meta.single_end
    def matefile_commands = meta.single_end
        ? ''
        : """
        mkdir ${prefix}_mate1 && mv ${prefix}_mate1.Chimeric.out.junction ${prefix}_mate1 && printf "${prefix}_mate1/${prefix}_mate1.Chimeric.out.junction" > mate1file
        mkdir ${prefix}_mate2 && mv ${prefix}_mate2.Chimeric.out.junction ${prefix}_mate2 && printf "${prefix}_mate2/${prefix}_mate2.Chimeric.out.junction" > mate2file
        """

    def mate_args = meta.single_end ? '' : '-mt1 @mate1file -mt2 @mate2file -Pi'

    """
    mkdir ${prefix} && mv ${prefix}.Chimeric.out.junction ${prefix} && printf "${prefix}/${prefix}.Chimeric.out.junction" > samplesheet
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
