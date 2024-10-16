process TARGETSCAN_DATABASE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(mature)

    output:
    tuple val(meta), path("mature.txt")  , emit: mature_txt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '1.3.4'
    """
    targetscan_format.sh $mature

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: $VERSION
    END_VERSIONS
    """
}
