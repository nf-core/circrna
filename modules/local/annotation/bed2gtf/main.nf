process ANNOTATION_BED2GTF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/785e4d9624ef7fb24cfee74f815d3e69540e99e0ae299cbd333e047e08706f7e/data' :
        'community.wave.seqera.io/library/polars_pyyaml:e53e9c9a38a99374' }"

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
