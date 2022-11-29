process FIND_CIRC_FILTER {
    tag "$meta.id"
    label "process_low"

    container 'barryd237/find_circ'

    input:
    tuple val(meta), path(bed)
    val bsj_reads

    output:
    tuple val(meta), path("${prefix}_find_circ_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_find_circ.bed")      , emit: matrix
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    grep circ ${prefix}.sites.bed | grep -v chrM | sum.py -2,3 | scorethresh.py -16 1 | scorethresh.py -15 2 | scorethresh.py -14 2 | scorethresh.py 7 ${bsj_reads} | scorethresh.py 8,9 35 | scorethresh.py -17 100000 >> ${prefix}.txt

    tail -n +2 ${prefix}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${prefix}_find_circ.bed

    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_find_circ.bed > ${prefix}_find_circ_circs.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find_circ: $VERSION
    END_VERSIONS
    """
}
