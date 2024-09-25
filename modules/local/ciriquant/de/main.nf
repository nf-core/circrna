process CIRIQUANT_DE {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/nicotru/ciriquant:1.0.4"

    input:
    tuple val(meta), path(library), path(expression), path(gene)

    output:
    tuple val(meta), path("${circ_path}"), emit: circ
    tuple val(meta), path("${gene_path}"), emit: gene
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    circ_path = "${prefix}.circ.csv"
    gene_path = "${prefix}.gene.csv"
    """
    CIRI_DE_replicate \\
        --lib ${library} \\
        --bsj ${expression} \\
        --gene ${gene} \\
        --out ${circ_path} \\
        --out2 ${gene_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
    END_VERSIONS
    """
}
