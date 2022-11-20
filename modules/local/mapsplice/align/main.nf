
process MAPSPLICE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mapsplice=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapsplice:2.2.1--py27h07887db_0':
        'quay.io/biocontainers/mapsplice:2.2.1--py27h07887db_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    path chromosomes
    path gtf

    output:
    tuple val(meta), path("${prefix}/fusions_raw.txt"), emit: raw_fusions
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 'v2.2.1'
    def gtf_prefix = gtf.toString() - ~/.gtf/
    def handleGzip_R1 = fastq[0].toString().endsWith('.gz') ? "gzip -d --force ${fastq[0]}" : ''
    def handleGzip_R2 = fastq[1].toString().endsWith('.gz') ? "gzip -d --force ${fastq[1]}" : ''
    def read1 = fastq[0].toString().endsWith('.gz') ? fastq[0].toString() - ~/.gz/ : fastq[0]
    def read2 = fastq[1].toString().endsWith('.gz') ? fastq[1].toString() - ~/.gz/ : fastq[1]
    """
    $handleGzip_R1
    $handleGzip_R2

    mapsplice.py \\
        -c $chromosomes \\
        -x $gtf_prefix \\
        -1 ${read1} \\
        -2 ${read2} \\
        -p ${task.cpus} \\
        --bam \\
        --gene-gtf $gtf \\
        -o $base \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapsplice: $VERSION
    END_VERSIONS
    """
}
