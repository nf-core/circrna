process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    tuple val(meta), path(gff)
    tuple val(meta2), path(genome)

    output:
    tuple val(meta), path("$outfile"), emit: output
    tuple val(meta), path("*.fasta") , emit: fasta, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def suffix        = task.ext.suffix ?: gff.extension
    def genome_string = genome ? "-g ${genome}" : ''
    outfile           = "${prefix}.${suffix}"
    """
    gffread \\
        ${genome_string} \\
        $args \\
        $gff > ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
