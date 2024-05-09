process SEGEMEHL_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(results)
    val(bsj_reads)

    output:
    tuple val(meta), path("${prefix}_segemehl_circs.bed"), emit: results
    tuple val(meta), path("${prefix}_segemehl.bed")      , emit: matrix
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.4'
    """
    grep ';C;' ${prefix}.sngl.bed | awk -v OFS="\\t" '{print \$1,\$2,\$3,\$6}' | sort | uniq -c | awk -v OFS="\\t" '{print \$2,\$3,\$4,\$5,\$1}' > ${prefix}_collapsed.bed

    awk -v OFS="\\t" -v BSJ=${bsj_reads} '{if(\$5>=BSJ) print \$0}' ${prefix}_collapsed.bed > ${prefix}_segemehl.bed

    awk -v OFS="\\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${prefix}_segemehl.bed > ${prefix}_segemehl_circs.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: $VERSION
        sort: \$(sort --version | head -n 1 | sed -e 's/sort (GNU coreutils) //g')
    END_VERSIONS
    """
}
