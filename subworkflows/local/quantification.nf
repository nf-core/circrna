include { PSIRC_QUANT                     } from './quantification_tools/psirc_quant'
include { CIRIQUANT                       } from './quantification_tools/ciriquant'
include { COMBINEBEDS_COUNTS as AGGREGATE } from '../../modules/local/combinebeds/counts'

workflow QUANTIFICATION {
    take:
    reads
    ch_gtf
    ch_fasta
    ch_transcriptome_fasta
    ch_transcriptome_gtf
    circ_annotation_bed
    circ_annotation_gtf
    ch_bsj_bed_per_sample_tool
    bootstrap_samples
    ch_phenotype
    ch_faidx
    bwa_index
    hisat2_index

    main:
    ch_versions = Channel.empty()
    ch_gene_counts = Channel.empty()
    ch_circ_counts = Channel.empty()
    ch_ciriquant = Channel.empty()
    ch_stringtie = Channel.empty()
    ch_rds = Channel.empty()

    tools_selected = params.quantification_tools.split(',').collect{it.trim().toLowerCase()}
    if (tools_selected.size() == 0) {
        error 'No tools selected for circRNA quantification.'
    }

    if (tools_selected.contains('psirc')) {
        PSIRC_QUANT(
            reads,
            ch_transcriptome_fasta,
            ch_transcriptome_gtf,
            circ_annotation_bed,
            circ_annotation_gtf,
            bootstrap_samples,
            ch_phenotype,
            ch_faidx
        )
        ch_gene_counts = ch_gene_counts
            .mix(PSIRC_QUANT.out.gene_counts.map{meta, counts -> [meta + [quantification: 'psirc'], counts]})
        ch_circ_counts = ch_circ_counts
            .mix(PSIRC_QUANT.out.circular_tx_counts.map{meta, counts -> [meta + [quantification: 'psirc'], counts]})
        ch_versions = ch_versions.mix(PSIRC_QUANT.out.versions)
        ch_rds = ch_rds.mix(PSIRC_QUANT.out.rds)
    }

    if (tools_selected.contains('ciriquant')) {
        CIRIQUANT(
            reads,
            circ_annotation_bed,
            ch_gtf,
            ch_fasta,
            bwa_index,
            hisat2_index
        )
        ch_versions = ch_versions.mix(CIRIQUANT.out.versions)
        ch_gene_counts = ch_gene_counts.mix(CIRIQUANT.out.gene_tpm)
        ch_circ_counts = ch_circ_counts.mix(CIRIQUANT.out.circ_cpm)
        ch_ciriquant = ch_ciriquant.mix(CIRIQUANT.out.raw)
        ch_stringtie = ch_stringtie.mix(CIRIQUANT.out.stringtie)
    }

    aggregations = ['sum', 'max']
    if (tools_selected.intersect(aggregations).size() > 0) {
        ch_aggregations = Channel.fromList(tools_selected.intersect(aggregations))

        AGGREGATE(
            ch_aggregations
                .map{agg -> [[id: agg], agg]}
                .combine(circ_annotation_bed.map{meta, bed -> bed})
                .combine(ch_bsj_bed_per_sample_tool.map{meta, bed -> bed}.collect().map{beds -> [beds]}),
            params.max_shift,
            params.consider_strand
        )
    }

    emit:
    gene      = ch_gene_counts
    circ      = ch_circ_counts
    ciriquant = ch_ciriquant
    stringtie = ch_stringtie
    rds       = ch_rds

    versions = ch_versions
}
