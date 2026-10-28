process CIRCTOOLS_RECONSTRUCT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f523ec9d03407da839bae40c76b9e49725d8b5273e9e167f7b5e081392f690d/data'
        : 'community.wave.seqera.io/library/circtools:2.0.5.post2--20ffd704d69916e2'}"

    input:
    tuple val(meta), path(bed), path(bam), path(bai), path(junction)
    tuple val(meta2), path(annotation)

    output:
    tuple val(meta), path("${prefix}.alternative_splicing.txt"), emit: alternative_splicing
    tuple val(meta), path("${prefix}.exon_counts.bed"), emit: exon_counts_bed
    tuple val(meta), path("${prefix}.exon_counts.txt"), emit: exon_counts_txt
    tuple val(meta), path("${prefix}.mate_status.txt"), emit: mate_status
    tuple val(meta), path("${prefix}.skipped_exons.bed"), emit: skipped_exons_bed
    tuple val(meta), path("${prefix}.skipped_exons.txt"), emit: skipped_exons_txt
    tuple val(meta), path("${prefix}/*.bam"), emit: bam, optional: true
    tuple val(meta), path("${prefix}.coverage_pictures/*.pdf"), emit: coverage_pictures, optional: true
    tuple val(meta), path("${prefix}.coverage_profiles/*.txt"), emit: coverage_profiles_txt, optional: true
    tuple val(meta), path("${prefix}.coverage_profiles/*.pdf"), emit: coverage_profiles_pdf, optional: true
    tuple val(meta), path("${prefix}.coverage_profiles/*.tsv"), emit: coverage_profiles_tsv, optional: true
    tuple val("${task.process}"), val('circtools'), eval('circtools -V'), emit: versions_circtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    circtools_package="\$(python3 -c 'import os, circtools; print(os.path.dirname(circtools.__file__))')"
    export PYTHONPATH="\$(dirname "\$circtools_package"):\$circtools_package:\$circtools_package/reconstruct:\${PYTHONPATH:-}"

    circtools reconstruct \
        --bamfile ${bam} \
        --annotation ${annotation} \
        --sampleName ${meta.id} \
        -D ${bed} \
        -J ${junction} \
        -T ./temp \
        -O . \
        -P ${task.cpus}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alternative_splicing.txt
    touch ${prefix}.exon_counts.bed
    touch ${prefix}.exon_counts.txt
    touch ${prefix}.mate_status.txt
    touch ${prefix}.skipped_exons.bed
    touch ${prefix}.skipped_exons.txt

    """
}
