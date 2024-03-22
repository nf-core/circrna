include { GNU_SORT as COMBINE_TRANSCRIPTOME_GTFS      } from '../../modules/nf-core/gnu/sort/main'
include { TRANSCRIPTOME                               } from '../../modules/local/quantification/transcriptome/main'
include { GAWK as MARK_CIRCULAR                       } from '../../modules/nf-core/gawk/main'
include { PSIRC_INDEX                                 } from '../../modules/local/psirc/index/main'
include { PSIRC_QUANT                                 } from '../../modules/local/psirc/quant/main'
include { CUSTOM_TX2GENE                              } from '../../modules/nf-core/custom/tx2gene/main'
include { TXIMETA_TXIMPORT                            } from '../../modules/nf-core/tximeta/tximport/main'
include { COMBINE_QUANTIFICATION as COMBINE_GENE_COUNTS    } from '../../modules/local/quantification/combine_quantification/main'
include { COMBINE_QUANTIFICATION as COMBINE_GENE_TPM       } from '../../modules/local/quantification/combine_quantification/main'
include { COMBINE_QUANTIFICATION as COMBINE_TX_COUNTS    } from '../../modules/local/quantification/combine_quantification/main'
include { COMBINE_QUANTIFICATION as COMBINE_TX_TPM       } from '../../modules/local/quantification/combine_quantification/main'

workflow QUANTIFICATION {
    take:
    ch_gtf
    ch_fasta
    counts
    reads
    circ_annotation_bed
    circ_annotation_gtf

    main:
    ch_versions = Channel.empty()

    COMBINE_TRANSCRIPTOME_GTFS(
        ch_gtf.mix(circ_annotation_gtf).map{meta, gtf -> gtf}.collect().map{[[id: "transcriptome"], it]},
    )

    TRANSCRIPTOME(COMBINE_TRANSCRIPTOME_GTFS.out.sorted, ch_fasta)
    MARK_CIRCULAR(TRANSCRIPTOME.out.transcriptome, [])

    ch_versions = ch_versions.mix(
        COMBINE_TRANSCRIPTOME_GTFS.out.versions,
        TRANSCRIPTOME.out.versions,
        MARK_CIRCULAR.out.versions
    )

    PSIRC_INDEX(MARK_CIRCULAR.out.output)
    PSIRC_QUANT(reads, PSIRC_INDEX.out.index.collect())

    CUSTOM_TX2GENE(
        COMBINE_TRANSCRIPTOME_GTFS.out.sorted,
        PSIRC_QUANT.out.directory.map{meta, quant -> quant}.collect().map{[[id: "quant"], it]},
        "kallisto",
        "gene_id",
        "gene_name"
    )

    TXIMETA_TXIMPORT(
        PSIRC_QUANT.out.directory,
        CUSTOM_TX2GENE.out.tx2gene,
        "kallisto"
    )

    ch_versions = ch_versions.mix(
        PSIRC_INDEX.out.versions,
        PSIRC_QUANT.out.versions,
        CUSTOM_TX2GENE.out.versions,
        TXIMETA_TXIMPORT.out.versions
    )

    ch_gene_counts = TXIMETA_TXIMPORT.out.counts_gene
    ch_gene_tpm = TXIMETA_TXIMPORT.out.tpm_gene

    COMBINE_GENE_COUNTS(
        ch_gene_counts.map{meta, counts -> counts}.collect().map{[[id: "gene_counts"], it]}
    )

    COMBINE_GENE_TPM(
        ch_gene_tpm.map{meta, tpm -> tpm}.collect().map{[[id: "gene_tpm"], it]}
    )

    COMBINE_TX_COUNTS(
        TXIMETA_TXIMPORT.out.counts_transcript.map{meta, counts -> counts}.collect().map{[[id: "tx_counts"], it]}
    )

    COMBINE_TX_TPM(
        TXIMETA_TXIMPORT.out.tpm_transcript.map{meta, tpm -> tpm}.collect().map{[[id: "tx_tpm"], it]}
    )

    ch_versions = ch_versions.mix(
        COMBINE_GENE_COUNTS.out.versions,
        COMBINE_GENE_TPM.out.versions
    )

    emit:
    gene_counts = COMBINE_GENE_COUNTS.out.combined
    gene_tpm    = COMBINE_GENE_TPM.out.combined

    versions        = ch_versions
}
