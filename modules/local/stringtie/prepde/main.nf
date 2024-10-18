process STRINGTIE_PREPDE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::stringtie=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    tuple val(meta), val(samples), path(gtfs)

    output:
    tuple val(meta), path("${prefix}_transcript_count_matrix.csv") , emit: transcript_matrix
    tuple val(meta), path("${prefix}_gene_count_matrix.csv")       , emit: gene_matrix

    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    samplesheet = [samples, gtfs]
        .transpose()
        .collect{ sample, gtf ->
            "${sample}\\t${gtf}" }.join('\\n')
    transcript_path = "${prefix}_transcript_count_matrix.csv"
    gene_path = "${prefix}_gene_count_matrix.csv"
    """
    echo -e "${samplesheet}" > samples.txt

    prepDE.py -i samples.txt \\
        -g ${gene_path} \\
        -t ${transcript_path} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
