process SJDB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sjdb)

    output:
    tuple val(meta), path("${prefix}.SJFile.tab"), emit: sjfile

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' !{sjdb} > ${prefix}.SJFile.tab
    """
}
