process STRINGTIE_PREPDE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/947b71729159919828a9d6560230b09524109152ed5edcc786e6f320ed368406/data'
        : 'community.wave.seqera.io/library/stringtie_prepde:bc1d541780932433'}"

    input:
    tuple val(meta), val(samples), path(gtfs)

    output:
    tuple val(meta), path("${prefix}_transcript_count_matrix.csv"), emit: transcript_matrix
    tuple val(meta), path("${prefix}_gene_count_matrix.csv"), emit: gene_matrix

    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version 2>&1'), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    samplesheet = [samples, gtfs]
        .transpose()
        .collect { sample, gtf ->
            "${sample}\\t${gtf}"
        }
        .join('\\n')
    transcript_path = "${prefix}_transcript_count_matrix.csv"
    gene_path = "${prefix}_gene_count_matrix.csv"
    """
    echo -e "${samplesheet}" > samples.txt

    prepDE.py -i samples.txt \\
        -g ${gene_path} \\
        -t ${transcript_path} \\
        ${args}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_transcript_count_matrix.csv
    touch ${prefix}_gene_count_matrix.csv

    """
}
