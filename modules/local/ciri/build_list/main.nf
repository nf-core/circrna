process BUILD_LIST {
    tag "${meta.id}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/polars_pyyaml:153427379c542734'
        : 'community.wave.seqera.io/library/polars_pyyaml:57e6c66da323f22b'}"

    input:
    tuple val(meta), path(ciri_annotation), path(consensus)
    val max_shift

    output:
    tuple val(meta), path("${prefix}.txt"), emit: list
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template "build_list.py"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
