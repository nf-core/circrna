process CIRIQUANT_DE {
    tag "$meta.id"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/nicotru/ciriquant:35228e054dcad0b0'
        : 'docker.io/nicotru/ciriquant@sha256:5e59c7a31355f1814a1703e9cb27d347b7bb04871606a24c1acdacfc02c627f8'}"

    input:
    tuple val(meta), path(library), path(expression), path(gene)

    output:
    tuple val(meta), path("${circ_path}"), emit: circ
    tuple val(meta), path("${gene_path}"), emit: gene
    tuple val("${task.process}"), val('ciriquant'), eval('CIRIquant --version 2>&1 | sed "s/CIRIquant //g"'), emit: versions_ciriquant, topic: versions
    tuple val("${task.process}"), val('r-base'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), emit: versions_r_base, topic: versions
    tuple val("${task.process}"), val('edger'), eval('Rscript -e "cat(as.character(packageVersion(\'edgeR\')))"'), emit: versions_edger, topic: versions

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

    if [ ! -f ${circ_path} ] || [ ! -f ${gene_path} ]; then
        echo "ERROR: CIRI_DE_replicate failed to produce ${circ_path} and ${gene_path}." >&2
        echo "edgeR most likely failed (e.g. too few samples or replicates per group); see .command.err and .command.log in the task directory." >&2
        exit 1
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    circ_path = "${prefix}.circ.csv"
    gene_path = "${prefix}.gene.csv"
    """
    touch ${circ_path}
    touch ${gene_path}

    """
}
