process JCCIRC_PREP {
    tag "${meta.id}"
    label 'process_low'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/polars_pyyaml:785c6487e524bef4'
        : 'community.wave.seqera.io/library/polars_pyyaml:92390185dab18544'}"

    input:
    tuple val(meta), path(bsj_annotation), path(bsj_reads)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: merged
    path "versions.yml", emit: versions
    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'prep.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        polars: \$(python -c "import polars as pl; print(pl.__version__)")
    END_VERSIONS
    """
}
