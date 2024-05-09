process DCC_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(txt)
    val bsj_reads

    output:
    tuple val(meta), path("${prefix}_dcc_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_dcc.bed")      , emit: matrix
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.4'
    """
    awk '{if(\$5 >= ${bsj_reads}) print \$0}' ${prefix}.txt > ${prefix}_dcc.filtered
    awk -v OFS="\\t" '{\$2-=1;print}' ${prefix}_dcc.filtered > ${prefix}_dcc.bed
    awk -v OFS="\\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_dcc.bed > ${prefix}_dcc_circs.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: $VERSION
    END_VERSIONS
    """
}
