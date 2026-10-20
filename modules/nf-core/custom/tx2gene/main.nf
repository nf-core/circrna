process CUSTOM_TX2GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5ce43c6bfb1120fe0a53ef6ca35b1a17597ef51ce071ae27cf5b68aa577c90b9/data'
        : 'community.wave.seqera.io/library/python:3.14.7--92bcd7a06268d464'}"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path ("quants/*")
    val quant_type
    val id
    val extra

    output:
    tuple val(meta), path("*tx2gene.tsv"), emit: tx2gene
    tuple val("${task.process}"), val('python'), eval('python --version | cut -d " " -f 2'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tx2gene.sh'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.tx2gene.tsv
    """
}
