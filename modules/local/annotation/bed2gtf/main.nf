process ANNOTATION_BED2GTF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_pyyaml:153427379c542734' :
        'community.wave.seqera.io/library/polars_pyyaml:57e6c66da323f22b' }"

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
