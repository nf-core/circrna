process CSVTK_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/673f8139d490e93cbfcc5e6f206ae2f516862c1b243a396cc34a0eef8eac8431/data'
        : 'community.wave.seqera.io/library/csvtk:0.31.0--07ef9a35e53e4e23'}"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("${prefix}.${out_extension}"), emit: csv
    tuple val("${task.process}"), val('csvtk'), eval('csvtk version | sed "s/csvtk v//"'), emit: versions_csvtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = args.contains('--out-delimiter "\t"') || args.contains('-D "\t"') || args.contains("-D \$'\t'") ? "tsv" : "csv"
    """
    csvtk \\
        summary \\
        $args \\
        --num-cpus $task.cpus \\
        --out-file ${prefix}.${out_extension} \\
        $csv

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    out_extension = args.contains('--out-delimiter "\t"') || args.contains('-D "\t"') || args.contains("-D \$'\t'") ? "tsv" : "csv"
    """
    touch ${prefix}.${out_extension}

    """
}
