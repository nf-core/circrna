process TXIMETA_TXIMETA {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/676a87512baa4aa94b5e3af327d07a75b170190cca324e25d4c9aa6d9b3f2fbd/data'
        : 'community.wave.seqera.io/library/tximeta_tximeta:91ac7b8609ffb2f2'}"

    input:
    tuple val(meta), path("quants/*")
    val quant_type

    output:
    tuple val(meta), path("${prefix}.rds"), emit: se
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template('tximeta.r')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: "\$(Rscript -e "cat(as.character(packageVersion('tximeta')))")"
    END_VERSIONS

    """
}
