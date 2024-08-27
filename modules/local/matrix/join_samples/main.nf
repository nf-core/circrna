process JOIN_SAMPLES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'oras://community.wave.seqera.io/library/polars_pyyaml:962a0cf7480258c7' :
            'community.wave.seqera.io/library/polars_pyyaml:ad93db0d7bcd508e' }"

    input:
    tuple val(meta), val(samples), path(matrices)

    output:
    tuple val(meta), path("${meta.id}.joined.tsv"), emit: joined
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    metacols = task.ext.metacols ?: "gene_id,gene_name"
    has_header = task.ext.has_header == null ? true : task.ext.has_header
    template 'join.py'

    stub:
    """
    touch ${meta.id}.joined.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
