process CIRCEXPLORER2_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(txt)
    val(bsj_reads)

    output:
    tuple val(meta), path("${prefix}_${meta.tool}_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_${meta.tool}.bed")      , emit: matrix
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.4'
    """
    awk '{if(\$13 >= ${bsj_reads}) print \$0}' ${prefix}.txt | awk -v OFS="\\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${prefix}_${meta.tool}.bed

    awk -v OFS="\\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_${meta.tool}.bed > ${prefix}_${meta.tool}_circs.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: $VERSION
    END_VERSIONS
    """
}
