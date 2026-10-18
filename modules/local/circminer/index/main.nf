process CIRCMINER_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/circminer:0.4.2--h5ca1c30_6'
        : 'quay.io/biocontainers/circminer:0.4.2--h5ca1c30_6'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('circminer'), emit: index
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir circminer
    # Do not let the indexer modify the staged FASTA symlink or its source.
    cp --dereference ${fasta} circminer/reference.fa
    (
        cd circminer
        circminer --index -r reference.fa -k 20 -t ${task.cpus} ${args}
    )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circminer: \$(circminer --version | sed -E 's/^CircMiner //')
    END_VERSIONS
    """

    stub:
    """
    mkdir circminer
    cp --dereference ${fasta} circminer/reference.fa
    touch circminer/reference.fa.packed.fa.index circminer/reference.fa.packed.fa.index.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circminer: 0.4.2
    END_VERSIONS
    """
}
