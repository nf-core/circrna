process CIRI_CIRILONG {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/openjdk:23.0.2--4046e6f926abce6e'
        : 'community.wave.seqera.io/library/openjdk:23.0.2--2fd1f5d679ee38ac'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    CIRI-long call -i ${reads} -r ${fasta} -a ${gtf} -p ${prefix} -t ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciri-long: \$(python3 -c "import CIRI_long; print(CIRI_long.__version__)")
    END_VERSIONS
    """
}
