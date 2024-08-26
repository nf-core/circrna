process CIRI_CIRIFULL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciri-full:2.1.2--hdfd78af_0' :
        'biocontainers/ciri-full:2.1.2--hdfd78af_0' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)
    
    output:
    tuple val(meta), path("${meta.id}.ciri"), emit: ciri

    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def VERSION = "2.1.2"
    """
    CIRI-full Pipeline -1 ${reads[0]} -2 ${reads[1]} -r ${fasta} -a ${gtf} -o ${meta.id} -d results -t ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri-full: ${VERSION}
    END_VERSIONS
    """
}