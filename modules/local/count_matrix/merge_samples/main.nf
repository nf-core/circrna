process MERGE_SAMPLES {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(beds)

    output:
    path("merged_counts.bed") , emit: counts_bed
    path("merged_counts.tsv") , emit: counts_tsv

    path "versions.yml"       , emit: versions

    script:
    template "merge_samples.py"
}
