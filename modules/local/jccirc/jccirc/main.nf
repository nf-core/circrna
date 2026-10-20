process JCCIRC_JCCIRC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d51b556a67811006e8517d14180695e4270c2b51b33b4c41682685b74a5916b2/data'
        : 'community.wave.seqera.io/library/jccirc:1.0.0--876376ce897b1bf5'}"

    input:
    tuple val(meta), path(reads), path(bsjs), path(denovo)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}.jccirc.log"), emit: log
    tuple val("${task.process}"), val('jccirc'), val('1.0.0'), emit: versions_jccirc, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    n_bsj=\$(wc -l ${bsjs} | awk '{print \$1}')
    max_allowed_cpus=${task.cpus}
    max_sensible_cpus=\$((n_bsj / 4))
    cpus=\$((max_sensible_cpus > max_allowed_cpus ? max_allowed_cpus : max_sensible_cpus))
    cpus=\$((cpus < 1 ? 1 : cpus))

    # JCcirc reports "... is empty!" and dies when no fragments are found. The exit code
    # of that die() depends on errno (e.g. 25 or 1), so the log is checked as well.
    # This is a valid empty result, not an error.
    status=0
    JCcirc -C ${bsjs} -G ${fasta} -F ${gtf} -P \$cpus --contig ${denovo} -O ${prefix} --read1 ${reads[0]} --read2 ${reads[1]} > ${prefix}.jccirc.log 2>&1 || status=\$?
    cat ${prefix}.jccirc.log

    if [ \$status -ne 0 ] && [ \$status -ne 25 ] && ! grep -q 'is empty!' ${prefix}.jccirc.log; then
        exit \$status
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.jccirc.log

    """
}
