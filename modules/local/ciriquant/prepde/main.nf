process CIRIQUANT_PREPDE {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/nicotru/ciriquant:1.0.4"

    input:
    tuple val(meta), val(samples), path(gtfs), val(conditions)

    output:
    tuple val(meta), path("${prefix}_library.tsv")   , emit: library
    tuple val(meta), path("${prefix}_annotation.tsv"), emit: annotation
    tuple val(meta), path("${prefix}_expression.tsv"), emit: expression
    tuple val(meta), path("${prefix}_ratio.tsv")     , emit: gene
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    samplesheet = [samples, gtfs, conditions]
        .transpose()
        .collect{ sample, gtf, condition ->
            "${sample}\\t${gtf}\\t${condition}" }.join('\\n')
    """
    echo -e "${samplesheet}" > samples.txt

    prep_CIRIquant -i samples.txt \\
        --lib ${prefix}_library.tsv \\
        --circ ${prefix}_annotation.tsv \\
        --bsj ${prefix}_expression.tsv \\
        --ratio ${prefix}_ratio.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
    END_VERSIONS
    """
}