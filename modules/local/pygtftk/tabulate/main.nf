process PYGTFTK_TABULATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'depot.galaxyproject.org/singularity/pygtftk:1.6.2--py39h4e691d4_2' :
        'biocontainers/pygtftk:1.6.2--py39h4e691d4_2' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("$outfile"), emit: output
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def suffix        = task.ext.suffix ?: gff.extension
    outfile           = "${prefix}.${suffix}"
    """
    gtftk tabulate \\
        $args \\
        -i $gtf > ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtftk: \$(gtftk -v | awk '{print substr(\$2, 2)}')
    END_VERSIONS
    """
}
