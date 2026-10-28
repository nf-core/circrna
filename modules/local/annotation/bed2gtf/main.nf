process ANNOTATION_BED2GTF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b584f817e89b571a59ff2e6c887863e005c4ba32057e188ef2868942dc9a9a4/data'
        : 'community.wave.seqera.io/library/python_polars:a7975bec8c9bd1e2'}"

    input:
    tuple val(meta), path(bed12), path(db_intersections)
    val exons_only

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template('bed2gtf.py')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: "\$(python3 -c 'import platform; print(platform.python_version())')"
        polars: "\$(python3 -c 'import polars; print(polars.__version__)')"
    END_VERSIONS
    """
}
