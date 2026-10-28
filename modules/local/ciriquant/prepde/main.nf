process CIRIQUANT_PREPDE {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://docker.io/nicotru/ciriquant:35228e054dcad0b0'
        : 'docker.io/nicotru/ciriquant@sha256:5e59c7a31355f1814a1703e9cb27d347b7bb04871606a24c1acdacfc02c627f8'}"

    input:
    tuple val(meta), val(samples), path(gtfs), val(conditions)

    output:
    tuple val(meta), path("${prefix}_library.tsv"), emit: library
    tuple val(meta), path("${prefix}_annotation.tsv"), emit: annotation
    tuple val(meta), path("${prefix}_expression.tsv"), emit: expression
    tuple val(meta), path("${prefix}_ratio.tsv"), emit: gene
    tuple val("${task.process}"), val('ciriquant'), eval('CIRIquant --version 2>&1 | sed "s/CIRIquant //g"'), emit: versions_ciriquant, topic: versions
    tuple val("${task.process}"), val('python'), eval('python --version | cut -d " " -f2'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    samplesheet = [samples, gtfs, conditions]
        .transpose()
        .collect { sample, gtf, condition ->
            "${sample}\\t${gtf}\\t${condition}"
        }
        .join('\\n')
    """
    echo -e "${samplesheet}" > samples.txt

    prep_CIRIquant -i samples.txt \\
        --lib ${prefix}_library.tsv \\
        --circ ${prefix}_annotation.tsv \\
        --bsj ${prefix}_expression.tsv \\
        --ratio ${prefix}_ratio.tsv \\
        ${args}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_library.tsv
    touch ${prefix}_annotation.tsv
    touch ${prefix}_expression.tsv
    touch ${prefix}_ratio.tsv

    """
}
