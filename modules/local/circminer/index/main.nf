process CIRCMINER_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4efcb00978baf72ee51c2f271993e85bec61a9fe1bcd67935632d2f44d941030/data'
        : 'community.wave.seqera.io/library/circminer:0.4.2--a7cbccf5a8b1afe7'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('circminer'), emit: index
    tuple val("${task.process}"), val('circminer'), eval('circminer --version | sed -E "s/^CircMiner //"'), emit: versions_circminer, topic: versions

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

    """

    stub:
    """
    mkdir circminer
    cp --dereference ${fasta} circminer/reference.fa
    touch circminer/reference.fa.packed.fa.index circminer/reference.fa.packed.fa.index.info

    """
}
