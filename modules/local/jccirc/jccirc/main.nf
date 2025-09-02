process JCCIRC_JCCIRC {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/jccirc:1.0.0--49e0f586e3d6ad30'
        : 'community.wave.seqera.io/library/jccirc:1.0.0--08c14a0bde3fbf0f'}"

    input:
    tuple val(meta), path(reads), path(bsjs), path(denovo)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.0.0"
    """
    n_bsj=\$(wc -l ${bsjs} | awk '{print \$1}')
    max_allowed_cpus=${task.cpus}
    max_sensible_cpus=\$((n_bsj / 4))
    cpus=\$((max_sensible_cpus > max_allowed_cpus ? max_allowed_cpus : max_sensible_cpus))

    JCcirc -C ${bsjs} -G ${fasta} -F ${gtf} -P \$cpus --contig ${denovo} -O ${prefix} --read1 ${reads[0]} --read2 ${reads[1]} || test \$? -eq 25

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jccirc: ${VERSION}
    END_VERSIONS
    """
}
