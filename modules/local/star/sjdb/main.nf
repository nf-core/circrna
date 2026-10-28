process SJDB {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dcdda74a12f8e69313d76c77a8612e16ffab53a8a42d979e62b1b3c4bede4c03/data'
        : 'community.wave.seqera.io/library/star_sjdb:517da11651f5309b'}"

    input:
    tuple val(meta), path(sjdb)
    val bsj_reads

    output:
    tuple val(meta), path("dataset.SJ.out.tab"), emit: sjtab
    tuple val("${task.process}"), val('mawk'), eval('mawk -W version 2>&1 | head -n 1'), emit: versions_mawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir tmp
    cat *.tab | awk -v BSJ=${bsj_reads} '(\$7 >= BSJ && \$6==0)' | cut -f1-6 | sort -T ./tmp/ | uniq > dataset.SJ.out.tab
    rm -rf tmp
    """

    stub:
    """
    touch dataset.SJ.out.tab

    """
}
