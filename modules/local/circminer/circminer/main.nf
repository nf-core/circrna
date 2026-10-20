process CIRCMINER_CIRCMINER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4efcb00978baf72ee51c2f271993e85bec61a9fe1bcd67935632d2f44d941030/data'
        : 'community.wave.seqera.io/library/circminer:0.4.2--a7cbccf5a8b1afe7'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}.circ_report"), emit: report
    tuple val("${task.process}"), val('circminer'), eval('circminer --version | sed -E "s/^CircMiner //"'), emit: versions_circminer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -- ${index}/*.packed.fa.index
    test "\$#" -eq 1 && test -f "\$1"
    reference=\${1%.packed.fa.index}
    test -f "\$reference"

    circminer -r "\$reference" -g ${gtf} -1 ${reads[0]} -2 ${reads[1]} -o ${prefix} -t ${task.cpus} ${args}
    # circMiner creates an empty report for a successful no-call run. Do not mask a crash
    # or another failure that did not create the report.
    test -f ${prefix}.circ_report

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.circ_report

    """
}
