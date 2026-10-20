process MAPSPLICE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/dba8f89114556c2bf1fe6bb94bef5253a754ca9854306facb1755ee221ce3b1c/data'
        : 'community.wave.seqera.io/library/mapsplice:046f5c757c7045ad'}"
    input:
    tuple val(meta), path(reads)
    path bowtie_index
    tuple val(meta2), path(chromosomes, stageAs: 'chromosomes/*')
    path gtf

    output:
    tuple val(meta), path("${prefix}/fusions_raw.txt"), emit: raw_fusions
    tuple val("${task.process}"), val('mapsplice'), val('2.2.1'), emit: versions_mapsplice, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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

        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/fusions_raw.txt

    """
}
