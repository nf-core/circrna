process SJDB {
    label 'process_single'

    conda "bioconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk%3A5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    path(sjdb)
    val(bsj_reads)

    output:
    path("dataset.SJ.out.tab"), emit: sjtab
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat *.tab | awk -v BSJ=${bsj_reads} '(\$7 >= BSJ && \$6==0)' | cut -f1-6 | sort | uniq > dataset.SJ.out.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
    END_VERSIONS
    """
}
