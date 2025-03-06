process ANNOTATION_BED2GTF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars:1.24.0--800cd3e4ff805434' :
        'community.wave.seqera.io/library/polars:1.24.0--2d2d323e8514e707' }"

    input:
    tuple val(meta), path(bed12), path(db_intersections)
    val exons_only

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: gtf

    path "versions.yml"                         , emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    suffix = task.ext.suffix ?: 'gtf'
    template 'bed2gtf.py'
}
