process INFER_STRAND {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::pysam=0.22.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.23.3--py311h5f69268_3'
        : 'quay.io/biocontainers/pysam:0.23.3--py311h5f69268_3'}"

    input:
    tuple val(meta), path(bed, stageAs: 'input/*')
    tuple val(refmeta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'infer_strand.py'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    cp ${bed} ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
