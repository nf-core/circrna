process FIND_CIRC {
    tag "${meta.id}"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebc79c1656071bb656528ff16510f3a9afb3dedffd14a08fd456a6c0ca19fd1f/data'
        : 'community.wave.seqera.io/library/find_circ:afd6068f3366fc58'}"

    input:
    tuple val(meta), path(anchors)
    tuple val(meta2), path(index)
    path fasta
    val strandedness

    output:
    tuple val(meta), path("${prefix}.sites.bed"), emit: bed
    tuple val(meta), path("${prefix}.sites.reads"), emit: reads
    tuple val(meta), path("${prefix}.sites.log"), emit: log
    tuple val("${task.process}"), val('bowtie2'), eval('bowtie2 --version 2>&1 | sed -n \'s/^bowtie2-align-s version //p\' | head -n 1'), emit: versions_bowtie2, topic: versions
    tuple val("${task.process}"), val('find_circ'), val('1.2'), emit: versions_find_circ, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ""
    args2 = task.ext.args2 ?: ""
    def strand_arg = strandedness == 'forward' || strandedness == 'reverse' ? "--stranded" : ""
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/.rev.1.bt2//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/.rev.1.bt2l//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        --threads ${task.cpus} \\
        --reorder \\
        --mm \\
        -D 20 \\
        --score-min=C,-15,0 \\
        -q \\
        -x \$INDEX \\
        ${args} \\
        -U ${anchors} | \\
        find_circ.py --genome=${fasta} ${strand_arg} ${args2} --prefix=${prefix} --stats=${prefix}.sites.log --reads=${prefix}.sites.reads > ${prefix}.sites.bed

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sites.bed
    touch ${prefix}.sites.reads
    touch ${prefix}.sites.log

    """
}
