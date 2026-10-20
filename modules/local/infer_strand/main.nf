process INFER_STRAND {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3db3d5f46165cb9c7e463f5b96a74b08612e2343dcc531f69f85c98cac2d95a0/data'
        : 'community.wave.seqera.io/library/infer_strand:84d996a7cde83840'}"

    input:
    tuple val(meta), path(bed, stageAs: 'input/*')
    tuple val(refmeta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('pysam'), eval('python -c "import pysam; print(pysam.__version__)"'), emit: versions_pysam, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template 'infer_strand.sh'

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    cp ${bed} ${prefix}.bed
    """
}
