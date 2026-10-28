process CIRI_CIRIVIS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a9c1f8c493dc731d72213644a92370b691ae4a6baf70eed5db0a5f6627428b42/data'
        : 'community.wave.seqera.io/library/openjdk:25.0.2--b6e848e158c9ff5b'}"

    input:
    tuple val(meta), path(anno), path(library_length), path(list)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}/${prefix}.list_circle.fa"), emit: fasta, optional: true
    tuple val(meta), path("${prefix}/${prefix}.list"), emit: list
    path "${prefix}/*.pdf", emit: pdf
    tuple val("${task.process}"), val('ciri-vis'), val('1.4.2'), emit: versions_ciri_vis, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def has_list = list ? true : false
    """
    if [ ${has_list} == false ]; then
        list=\$(for file in ${anno}; do
            tail -n +2 \$file | cut -f1 | sort -u
        done)

        LIST_FILE="circle_ids.list"
        printf "%s\\n" \$list > \$LIST_FILE
    else
        LIST_FILE=${list}
    fi

    java -jar \$(command -v CIRI_vis_v1.4.2.jar) -i ${anno} -l ${library_length} -r ${fasta} -d ${prefix} -o ${prefix} -list \$LIST_FILE ${args}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.list_circle.fa
    touch ${prefix}/${prefix}.list
    touch ${prefix}/${prefix}.pdf

    """
}
