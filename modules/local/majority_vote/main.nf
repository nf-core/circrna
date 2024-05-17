process MAJORITY_VOTE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
		'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(miranda_data)
	path(targetscan_data)

    output:
    tuple val(meta), path("${meta.id}.majority.tsv.gz")               , emit: tsv 
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'majority.R'
}
