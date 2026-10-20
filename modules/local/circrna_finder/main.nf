process CIRCRNA_FINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/ee41e12abc2b47b1ce60ab8b6448961c701c00582098810066b5bae2c206bf03/data'
        : 'community.wave.seqera.io/library/circrna_finder:b96e4b698e6cee7c'}"

    input:
    tuple val(meta), path(star_input, stageAs: 'input/')

    output:
    tuple val(meta), path("${prefix}.filteredJunctions.bed"), emit: results
    tuple val("${task.process}"), val('awk'), eval('awk --version | head -n1 | cut -d\' \' -f3 | sed \'s/,//g\''), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('cat'), eval('cat --version | head -n1 | sed \'s/cat (GNU coreutils) //\''), emit: versions_cat, topic: versions
    tuple val("${task.process}"), val('circRNA_finder'), val('1.2'), emit: versions_circrna_finder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    postProcessStarAlignment.pl ${args} --starDir input/ --outDir ./

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filteredJunctions.bed

    """
}
