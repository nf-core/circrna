process ANNOTATION_BED2GTF {
    tag "$meta.id"
    label 'process_medium'
    shell '/bin/bash'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b584f817e89b571a59ff2e6c887863e005c4ba32057e188ef2868942dc9a9a4/data'
        : 'community.wave.seqera.io/library/python_polars:a7975bec8c9bd1e2'}"

    input:
    tuple val(meta), path(bed12), path(db_intersections)
    val exons_only

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('polars'), eval("python -c 'import polars; print(polars.__version__)'"), emit: versions_polars, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'bed2gtf.sh'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.gtf

    """
}
