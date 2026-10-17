process COMBINEBEDS_COUNTS {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/785e4d9624ef7fb24cfee74f815d3e69540e99e0ae299cbd333e047e08706f7e/data' :
        'community.wave.seqera.io/library/polars_pyyaml:e53e9c9a38a99374' }"

    input:
    tuple val(meta), val(aggregation), path(candidates), path(beds)
    val(max_shift)
    val(consider_strand)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: combined, optional: true
    path "versions.yml"                   , emit: versions

    script:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    template "counts.py"

    stub:
    prefix      = task.ext.prefix      ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
