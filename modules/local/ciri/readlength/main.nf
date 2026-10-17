process CIRI_READLENGTH {
    tag "${meta.id}"

    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/polars_pyyaml:785c6487e524bef4'
        : 'community.wave.seqera.io/library/polars_pyyaml:92390185dab18544'}"

    input:
    tuple val(meta), path(lengths)

    output:
    tuple val(meta), path("*.txt*"), emit: length
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    template('readlength.py')

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    echo 100 > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
