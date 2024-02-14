process SJDB {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(sjdb)
    val(bsj_reads)

    output:
    tuple val(meta), path("dataset.SJ.out.tab"), emit: sjtab
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '1.3.4'
    """
    mkdir tmp
    cat *.tab | awk -v BSJ=${bsj_reads} '(\$7 >= BSJ && \$6==0)' | cut -f1-6 | sort -T ./tmp/ | uniq > dataset.SJ.out.tab
    rm -rf tmp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: $VERSION
    END_VERSIONS
    """
}
