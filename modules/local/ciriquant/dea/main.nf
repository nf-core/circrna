process CIRIQUANT_DEA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ciriquant=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.2--pyhdfd78af_2' :
        'biocontainers/ciriquant:1.1.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), val(samples), path(gtfs), val(conditions)

    output:
    tuple val(meta), path("${library_path}"), path("${annotation_path}"), path("${expression_path}"), path("${ratio_path}"), emit: results
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    samplesheet = [samples, gtfs, conditions]
        .transpose()
        .collect{ sample, gtf, condition ->
            "${sample}\t${gtf}\t${condition}" }.join('\n')
    samplesheet_path = "${prefix}_samplesheet.tsv"
    library_path = "${prefix}_library.tsv"
    annotation_path = "${prefix}_annotation.tsv"
    expression_path = "${prefix}_expression.tsv"
    ratio_path = "${prefix}_ratio.tsv"
    """
    echo -e "${samplesheet}" > ${samplesheet_path}

    prep_CIRIquant -i ${samplesheet_path} \\
        --lib ${library_path} \\
        --circ ${annotation_path} \\
        --bsj ${expression_path} \\
        --ratio ${ratio_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
    END_VERSIONS
    """
}