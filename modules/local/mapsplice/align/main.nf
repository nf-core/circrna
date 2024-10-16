process MAPSPLICE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mapsplice=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapsplice:2.2.1--py27h07887db_0':
        'biocontainers/mapsplice:2.2.1--py27h07887db_0' }"

    input:
    tuple val(meta), path(reads)
    path bowtie_index
    tuple val(meta2), path(chromosomes, stageAs: 'chromosomes/*')
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
    if(meta.single_end){
        def handleGzip_R1 = reads[0].getExtension() == 'gz' ? "gzip -d -f ${reads[0]}" : ''
        def read1 = reads[0].getExtension() == 'gz' ? reads[0].toString() - ~/.gz/ : reads[0]
        """
        $handleGzip_R1

        mapsplice.py \\
            -c chromosomes \\
            -x $gtf_prefix \\
            -1 ${read1} \\
            -p ${task.cpus} \\
            --bam \\
            --gene-gtf $gtf \\
            -o $prefix \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mapsplice: $VERSION
        END_VERSIONS
        """
    } else {
        def handleGzip_R1 = reads[0].getExtension() == 'gz' ? "gzip -d -f ${reads[0]}" : ''
        def handleGzip_R2 = reads[1].getExtension() == 'gz' ? "gzip -d -f ${reads[1]}" : ''
        def read1 = reads[0].getExtension() == 'gz' ? reads[0].toString() - ~/.gz/ : reads[0]
        def read2 = reads[1].getExtension() == 'gz' ? reads[1].toString() - ~/.gz/ : reads[1]
        """
        $handleGzip_R1
        $handleGzip_R2

        mapsplice.py \\
            -c chromosomes \\
            -x $gtf_prefix \\
            -1 ${read1} \\
            -2 ${read2} \\
            -p ${task.cpus} \\
            --bam \\
            --gene-gtf $gtf \\
            -o $prefix \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mapsplice: $VERSION
        END_VERSIONS
        """
    }
}
