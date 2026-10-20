process PYGTFTK_TABULATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e3b69f00edfd9e23a400689d62645ed0588d0f0504cc589536df8bab2ded673/data'
        : 'community.wave.seqera.io/library/pygtftk_tabulate:db1d4cb61817f9b1'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("$outfile"), emit: table
    tuple val("${task.process}"), val('gtftk'), val('1.6.2'), emit: versions_gtftk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile    = "${prefix}.circ.tsv"
    """
    export HOME="\$PWD"

    gtftk tabulate \\
        $args \\
        -i $gtf | \\
        grep -v '^tabulate()' > ${outfile}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    outfile    = "${prefix}.circ.tsv"
    """
    export HOME="\$PWD"

    touch ${outfile}

    """
}
