process DCC_FILTER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(txt)
    val bsj_reads

    output:
    tuple val(meta), path("${prefix}_dcc_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_dcc.bed")      , emit: matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '{if(\$5 >= ${bsj_reads}) print \$0}' ${prefix}.txt > ${prefix}_dcc.filtered
    awk -v OFS="\t" '{\$2-=1;print}' ${prefix}_dcc.filtered > ${prefix}_dcc.bed
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_dcc.bed > ${prefix}_dcc_circs.bed
    """
}
