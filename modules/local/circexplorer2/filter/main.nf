process CIRCEXPLORER2_FILTER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(txt)
    val(bsj_reads)

    output:
    tuple val(meta), path("${prefix}_circexplorer2_circs.bed"), emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{if(\$13 >= ${bsj_reads}) print \$0}' ${prefix}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${prefix}_circexplorer2.bed

    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_circexplorer2.bed > ${prefix}_circexplorer2_circs.bed
    """
}
