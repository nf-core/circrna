process CIRI_CIRI2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/ciri-full_samtools:858c656b1f4e9732' :
        'community.wave.seqera.io/library/ciri-full_samtools:754497818ad973a6' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)
    
    output:
    tuple val(meta), path("${meta.id}.ciri"), emit: ciri

    path "versions.yml", emit: versions
    
    script:
    def sam = "${bam.baseName}.sam"
    def args = task.ext.args ?: ''
    """
    samtools view -h ${bam} -o ${sam}
    CIRI -I ${sam} -O ${meta.id}.ciri -F ${fasta} -A ${gtf} --log ${meta.id}.log -T ${task.cpus} ${args}
    rm ${sam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri: \$(echo \$(CIRI --help 2>&1) | sed 's/^.*Version: //; s/ Contact.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}