process CIRIQUANT_FILTER {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(gtf)
    val bsj_reads

    output:
    tuple val(meta), path("${prefix}_ciriquant_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_ciriquant.bed")      , emit: matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -v "#" ${prefix}.gtf | awk '{print \$14}' | cut -d '.' -f1 > counts
    grep -v "#" ${prefix}.gtf | awk -v OFS="\t" '{print \$1,\$4,\$5,\$7}' > ${prefix}.tmp
    paste ${prefix}.tmp counts > ${prefix}_unfilt.bed

    awk '{if(\$5 >= ${bsj_reads}) print \$0}' ${prefix}_unfilt.bed > ${prefix}_filt.bed
    grep -v '^\$' ${prefix}_filt.bed > ${prefix}_ciriquant

    awk -v OFS="\t" '{\$2-=1;print}' ${prefix}_ciriquant > ${prefix}_ciriquant.bed
    rm ${prefix}.gtf

    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_ciriquant.bed > ${prefix}_ciriquant_circs.bed
    """
}
