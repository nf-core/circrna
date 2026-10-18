process CIRI3_UC_LIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/polars_pyyaml:ad93db0d7bcd508e'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${prefix}.ciri3.list"), emit: list
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'uc_list.py'
    """
    python3 uc_list.py ${bed} ${prefix}.ciri3.list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ciri3.list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
