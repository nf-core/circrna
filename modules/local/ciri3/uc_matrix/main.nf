process CIRI3_UC_MATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/polars_pyyaml:ad93db0d7bcd508e'

    input:
    tuple val(meta), path(circ_list)
    tuple val(meta2), path(bsj_matrix), path(fsj_matrix)

    output:
    tuple val(meta2), path("${prefix}.bsj.tsv"), emit: bsj
    tuple val(meta2), path("${prefix}.fsj.tsv"), emit: fsj
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta2.id}"
    template 'uc_matrix.py'
    """
    python3 uc_matrix.py ${circ_list} ${bsj_matrix} ${fsj_matrix} ${prefix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta2.id}"
    """
    printf 'circRNA_ID\\tcount\\n' > ${prefix}.bsj.tsv
    printf 'circRNA_ID\\tcount\\n' > ${prefix}.fsj.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
