process JCCIRC_JCCIRC {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/jccirc:1.0.0--49e0f586e3d6ad30'
        : 'community.wave.seqera.io/library/jccirc:1.0.0--08c14a0bde3fbf0f'}"

    input:
    tuple val(meta), path(reads), path(bsjs), path(denovo)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}_ro1_align.txt"), emit: align
    tuple val(meta), path("${prefix}_ro1.fq.gz"), emit: fastq
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "2.1.2"
    """
    JCcirc -C ${bsjs} -G ${fasta} -F ${gtf} -P ${task.cpus} --contig ${denovo} -O ${prefix} --read1 ${reads[0]} --read2 ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cirifull: ${VERSION}
    END_VERSIONS
    """
}
