include { GAWK as MARK_CIRCULAR                  } from '../../../modules/nf-core/gawk'
include { PSIRC_INDEX                            } from '../../../modules/local/psirc/index'
include { PSIRC_QUANT as RUN_PSIRC_QUANT         } from '../../../modules/local/psirc/quant'
include { CUSTOM_TX2GENE                         } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT                       } from '../../../modules/nf-core/tximeta/tximport'
include { TXIMETA_TXIMETA                        } from '../../../modules/local/tximeta/tximeta'
include { MERGE_EXPERIMENTS                      } from '../../../modules/local/quantification/merge_experiments'
include { CSVTK_JOIN as JOIN_GENE_COUNTS         } from '../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_GENE_TPM            } from '../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_TX_COUNTS           } from '../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_TX_TPM              } from '../../../modules/nf-core/csvtk/join'
include { SPLIT_TYPES as SPLIT_TYPES_COUNTS      } from '../../../modules/local/quantification/split_types'
include { SPLIT_TYPES as SPLIT_TYPES_TPM         } from '../../../modules/local/quantification/split_types'

workflow PSIRC_QUANT {
    take:
    reads
    ch_transcriptome_fasta
    ch_transcriptome_gtf
    circ_annotation_bed
    circ_annotation_gtf
    bootstrap_samples
    ch_phenotype
    ch_faidx

    main:
    ch_versions = Channel.empty()

    MARK_CIRCULAR(ch_transcriptome_fasta, [])
    ch_versions = ch_versions.mix(MARK_CIRCULAR.out.versions)

    PSIRC_INDEX(MARK_CIRCULAR.out.output)
    RUN_PSIRC_QUANT(reads, PSIRC_INDEX.out.index.collect(), MARK_CIRCULAR.out.output, ch_faidx, bootstrap_samples)

    CUSTOM_TX2GENE(
        ch_transcriptome_gtf,
        RUN_PSIRC_QUANT.out.directory.map{meta, quant -> quant}.collect().map{[[id: "quant"], it]},
        "kallisto",
        "gene_id",
        "gene_name"
    )

    TXIMETA_TXIMETA(
        RUN_PSIRC_QUANT.out.directory,
        "kallisto"
    )

    TXIMETA_TXIMPORT(
        RUN_PSIRC_QUANT.out.directory,
        CUSTOM_TX2GENE.out.tx2gene,
        "kallisto"
    )

    ch_versions = ch_versions.mix(
        PSIRC_INDEX.out.versions,
        RUN_PSIRC_QUANT.out.versions,
        CUSTOM_TX2GENE.out.versions,
        TXIMETA_TXIMETA.out.versions,
        TXIMETA_TXIMPORT.out.versions
    )

    JOIN_GENE_COUNTS(
        TXIMETA_TXIMPORT.out.counts_gene.map{meta, counts -> counts}.collect().map{[[id: "gene_counts"], it]}
    )

    JOIN_GENE_TPM(
        TXIMETA_TXIMPORT.out.tpm_gene.map{meta, tpm -> tpm}.collect().map{[[id: "gene_tpm"], it]}
    )

    JOIN_TX_COUNTS(
        TXIMETA_TXIMPORT.out.counts_transcript.map{meta, counts -> counts}.collect().map{[[id: "tx_counts"], it]}
    )

    JOIN_TX_TPM(
        TXIMETA_TXIMPORT.out.tpm_transcript.map{meta, tpm -> tpm}.collect().map{[[id: "tx_tpm"], it]}
    )

    SPLIT_TYPES_COUNTS(
        JOIN_TX_COUNTS.out.csv
    )

    SPLIT_TYPES_TPM(
        JOIN_TX_TPM.out.csv
    )


    MERGE_EXPERIMENTS(
        TXIMETA_TXIMETA.out.se.map{meta, se -> se}.collect().map{[[id: "experiments"], it]},
        ch_phenotype.ifEmpty([[], []]),
        ch_transcriptome_gtf,
        JOIN_TX_TPM.out.csv
    )

    ch_versions = ch_versions.mix(
        JOIN_GENE_COUNTS.out.versions,
        JOIN_GENE_TPM.out.versions,
        JOIN_TX_COUNTS.out.versions,
        JOIN_TX_TPM.out.versions,
        SPLIT_TYPES_COUNTS.out.versions,
        SPLIT_TYPES_TPM.out.versions,
        MERGE_EXPERIMENTS.out.versions
    )

    emit:
    se                 = MERGE_EXPERIMENTS.out.merged
    rds                = MERGE_EXPERIMENTS.out.merged
    gene_counts        = JOIN_GENE_COUNTS.out.csv
    gene_tpm           = JOIN_GENE_TPM.out.csv
    tx_counts          = JOIN_TX_COUNTS.out.csv
    tx_tpm             = JOIN_TX_TPM.out.csv
    linear_tx_counts   = SPLIT_TYPES_COUNTS.out.linear
    linear_tx_tpm      = SPLIT_TYPES_TPM.out.linear
    circular_tx_counts = SPLIT_TYPES_COUNTS.out.circular
    circular_tx_tpm    = SPLIT_TYPES_TPM.out.circular

    versions           = ch_versions
}