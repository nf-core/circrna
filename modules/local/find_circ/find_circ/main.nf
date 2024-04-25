process FIND_CIRC {
    tag "$meta.id"
    label "process_high"

    conda "bioconda::find_circ=1.2 bioconda::bowtie2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c27e472038a09e49d9147bc52903e12836302c12:60ffb3b15a2c40c669f8d38382b1e6e4b065f5e4-0' :
        'biocontainers/mulled-v2-c27e472038a09e49d9147bc52903e12836302c12:60ffb3b15a2c40c669f8d38382b1e6e4b065f5e4-0' }"

    input:
    tuple val(meta), path(anchors)
    tuple val(meta2), path(index)
    path fasta

    output:
    tuple val(meta), path("${prefix}.sites.bed"), emit: bed
    path("${prefix}.sites.reads")               , emit: reads
    path("${prefix}.sites.log")                 , emit: log
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2'
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
        find_circ.py  --genome=$fasta --prefix=${prefix} --stats=${prefix}.sites.log --reads=${prefix}.sites.reads > ${prefix}.sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        find_circ: $VERSION
    END_VERSIONS
    """
}
