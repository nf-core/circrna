process TARGETSCAN {
    tag "$meta.id"
    label 'process_medium'

    container 'barryd237/targetscan'

    input:
    tuple val(meta), path(fasta)
    path(mature_txt)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "7.0"
    """
    ##format for targetscan
    cat $fasta | grep ">" | sed 's/>//g' > id
    cat $fasta | grep -v ">" > seq
    paste id seq | awk -v OFS="\t" '{print \$1, "0000", \$2}' > ${prefix}_ts.txt
    # run targetscan
    targetscan_70.pl mature.txt ${prefix}_ts.txt ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        targetscan: $VERSION
    END_VERSIONS
    """
}
