process CIRIQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ciriquant=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.2--pyhdfd78af_2' :
        'quay.io/biocontainers/ciriquant:1.1.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    path yml

    output:
    tuple val(meta), path("${prefix}/${prefix}.gtf"), emit: gtf
    path("${prefix}")                               , emit: results
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.0'
    """
    CIRIquant \\
        -t ${task.cpus} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --config $yml \\
        --no-gene \\
        -o ${prefix} \\
        -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        ciriquant : \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        stringtie: \$(stringtie --version 2>&1)
        hisat2: $VERSION
    END_VERSIONS
    """
}
