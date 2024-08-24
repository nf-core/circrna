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
    tuple val(meta), path("${samplesheet_path}"), emit: samplesheet
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
    """
    echo -e "${samplesheet}" > ${samplesheet_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
    END_VERSIONS
    """
}