process CIRCMINER_CIRCMINER {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/circminer:0.4.2--h5ca1c30_6'
        : 'quay.io/biocontainers/circminer:0.4.2--h5ca1c30_6'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}.circ_report"), emit: report
    path 'versions.yml', emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circminer: \$(circminer --version | sed -E 's/^CircMiner //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.circ_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circminer: 0.4.2
    END_VERSIONS
    """
}
