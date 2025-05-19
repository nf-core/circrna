process CIRI_CIRIVIS {
    tag "${meta.id}"
    label 'process_low'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/openjdk:23.0.2--4046e6f926abce6e'
        : 'community.wave.seqera.io/library/openjdk:23.0.2--2fd1f5d679ee38ac'}"

    input:
    tuple val(meta), path(anno), path(library_length), path(list)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}/${prefix}.list_circle.fa"), emit: fasta, optional: true
    tuple val(meta), path("${prefix}/${prefix}.list"), emit: list
    path "${prefix}/*.pdf", emit: pdf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.4.2"
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

    java -jar ${moduleDir}/CIRI_vis_v${VERSION}.jar -i ${anno} -l ${library_length} -r ${fasta} -d ${prefix} -o ${prefix} -list \$LIST_FILE ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri-vis: ${VERSION}
    END_VERSIONS
    """
}
