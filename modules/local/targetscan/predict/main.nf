process TARGETSCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/2854bd21d12ba23d6df21cbda53ec46f5c5b860720ef48bf4ae3b16c8a97e92b/data'
        : 'community.wave.seqera.io/library/targetscan:9c82f6ff6777613b'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(mature_txt)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('awk'), eval('awk --version | head -n1 | cut -d\' \' -f3 | sed \'s/,//g\''), emit: versions_awk, topic: versions
    tuple val("${task.process}"), val('targetscan'), val('7.0'), emit: versions_targetscan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ##format for targetscan
    cat ${fasta} | grep ">" | sed 's/>//g' > id
    cat ${fasta} | grep -v ">" > seq
    paste id seq | awk -v OFS="\\t" '{print \$1, "0000", \$2}' > ${prefix}_ts.txt
    # run targetscan
    targetscan_70.pl ${mature_txt} ${prefix}_ts.txt ${prefix}.txt

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    """
}
