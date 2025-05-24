process CIRCTOOLS_RECONSTRUCT {
    tag "${meta.id}"
    label 'process_high'

    conda "environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/circtools:2.0--b33a23b50d9f0697'
        : 'community.wave.seqera.io/library/circtools:2.0--f5bc60d7f93fefae'}"

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
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    circtools reconstruct \
        --bamfile ${bam} \
        --annotation ${annotation} \
        --sampleName ${meta.id} \
        -D ${bed} \
        -J ${junction} \
        -T ./temp \
        -O . \
        -P ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circtools: \$(circtools -V)
    END_VERSIONS
    """
}
