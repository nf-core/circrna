process FIND_CIRC {
    tag "$meta.id"
    label "process_high"

    container 'barryd237/find_circ'

    input:
    tuple val(meta), path(anchors)
    tuple val(meta2), path(index)
    path fasta
    path chromosomes

    output:
    tuple val(meta), path("${prefix}.sites.bed"), emit: bed
    path("${prefix}.sites.reads")               , emit: reads
    path("${prefix}.sites.log")                 , emit: log
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/.rev.1.bt2//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/.rev.1.bt2l//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        --threads $task.cpus \\
        --reorder \\
        --mm \\
        -D 20 \\
        --score-min=C,-15,0 \\
        -q \\
        -x \$INDEX \\
        -U $anchors | \\
        find_circ.py  -G $chromosomes -p ${prefix} -s ${prefix}.sites.log > ${prefix}.sites.bed 2> ${prefix}.sites.reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        find_circ: $VERSION
    END_VERSIONS
    """
}
