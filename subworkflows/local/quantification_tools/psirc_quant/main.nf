include { GAWK as MARK_CIRCULAR             } from '../../../../modules/nf-core/gawk'
include { PSIRCQUANT_INDEX as INDEX         } from '../../../../modules/local/psircquant/index'
include { PSIRCQUANT_QUANT as QUANT         } from '../../../../modules/local/psircquant/quant'
include { CUSTOM_TX2GENE                    } from '../../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT                  } from '../../../../modules/nf-core/tximeta/tximport'
include { TXIMETA_TXIMETA                   } from '../../../../modules/local/tximeta/tximeta'
include { MERGE_EXPERIMENTS                 } from '../../../../modules/local/quantification/merge_experiments'
include { CSVTK_JOIN as JOIN_GENE_COUNTS    } from '../../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_GENE_TPM       } from '../../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_TX_COUNTS      } from '../../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_TX_TPM         } from '../../../../modules/nf-core/csvtk/join'
include { SPLIT_TYPES as SPLIT_TYPES_COUNTS } from '../../../../modules/local/quantification/split_types'
include { SPLIT_TYPES as SPLIT_TYPES_TPM    } from '../../../../modules/local/quantification/split_types'

workflow PSIRC_QUANT {
    take:
    reads
    ch_transcriptome_fasta
    ch_transcriptome_gtf
    bootstrap_samples
    ch_phenotype
    ch_faidx

    main:

    MARK_CIRCULAR(ch_transcriptome_fasta.map { meta, fasta -> [meta, fasta, 'marked.fasta'] }, [], false)

    INDEX(MARK_CIRCULAR.out.output)

    QUANT(reads, INDEX.out.index.collect(), ch_transcriptome_gtf, ch_faidx, bootstrap_samples)

    CUSTOM_TX2GENE(
        ch_transcriptome_gtf,
        QUANT.out.directory.map { _meta, quant -> quant }.collect().map { quants -> [[id: "quant"], quants] },
        "kallisto",
        "gene_id",
        "gene_name",
    )

    TXIMETA_TXIMETA(
        QUANT.out.directory,
        "kallisto",
    )

    TXIMETA_TXIMPORT(
        QUANT.out.directory,
        CUSTOM_TX2GENE.out.tx2gene,
        "kallisto",
    )

    JOIN_GENE_COUNTS(
        TXIMETA_TXIMPORT.out.counts_gene.map { _meta, file -> file }.collect(sort: { file -> file.name }).map { files -> [[id: "gene_counts"], files] }
    )

    JOIN_GENE_TPM(
        TXIMETA_TXIMPORT.out.tpm_gene.map { _meta, file -> file }.collect(sort: { file -> file.name }).map { files -> [[id: "gene_tpm"], files] }
    )

    JOIN_TX_COUNTS(
        TXIMETA_TXIMPORT.out.counts_transcript.map { _meta, file -> file }.collect(sort: { file -> file.name }).map { files -> [[id: "tx_counts"], files] }
    )

    JOIN_TX_TPM(
        TXIMETA_TXIMPORT.out.tpm_transcript.map { _meta, file -> file }.collect(sort: { file -> file.name }).map { files -> [[id: "tx_tpm"], files] }
    )

    SPLIT_TYPES_COUNTS(
        JOIN_TX_COUNTS.out.out_file
    )

    SPLIT_TYPES_TPM(
        JOIN_TX_TPM.out.out_file
    )

    MERGE_EXPERIMENTS(
        TXIMETA_TXIMETA.out.se.map { _meta, se -> se }.collect(sort: { se -> se.name }).map { experiments -> [[id: "experiments"], experiments] },
        ch_phenotype.ifEmpty([[], []]),
        ch_transcriptome_gtf,
        JOIN_TX_TPM.out.out_file,
    )

    emit:
    se                 = MERGE_EXPERIMENTS.out.merged
    rds                = MERGE_EXPERIMENTS.out.merged
    gene_counts        = JOIN_GENE_COUNTS.out.out_file
    gene_tpm           = JOIN_GENE_TPM.out.out_file
    tx_counts          = JOIN_TX_COUNTS.out.out_file
    tx_tpm             = JOIN_TX_TPM.out.out_file
    linear_tx_counts   = SPLIT_TYPES_COUNTS.out.linear
    linear_tx_tpm      = SPLIT_TYPES_TPM.out.linear
    circular_tx_counts = SPLIT_TYPES_COUNTS.out.circular
    circular_tx_tpm    = SPLIT_TYPES_TPM.out.circular
}
