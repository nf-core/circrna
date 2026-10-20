process PSIRCQUANT_QUANT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e5afa9a01488ef93da1c7d72dfb0f0c9d885e5933aa67aa4ccb52f3c98eeb0c/data'
        : 'community.wave.seqera.io/library/psirc:1.0.0--4c24a8f727283097'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    tuple val(meta4), path(chrom_sizes)
    val bootstrap_samples

    output:
    tuple val(meta), path("${meta.id}"), emit: directory
    tuple val("${task.process}"), val('psirc-quant'), eval('psirc-quant version | sed -n "s/^psirc-quant, version \\([0-9.]*\\).*$/\\1/p"'), emit: versions_psircquant, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def single_end = meta.single_end ? "--single -l 76 -s 20" : ""
    def genomebam = gtf ? "--genomebam -g ${gtf}" : ""
    def chromosomes = chrom_sizes ? "-c ${chrom_sizes}" : ""
    """
    psirc-quant quant -t ${task.cpus} -i ${index} -o ${meta.id} ${single_end} ${reads} -b ${bootstrap_samples} ${genomebam} ${chromosomes}

    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/abundance.h5
    touch ${meta.id}/abundance.tsv
    touch ${meta.id}/run_info.json

    """
}
