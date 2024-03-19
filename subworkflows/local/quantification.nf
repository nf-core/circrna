include { TRANSCRIPTOME as LINEAR_TRANSCRIPTOME       } from '../../modules/local/quantification/transcriptome/main'
include { BEDTOOLS_GETFASTA as CIRCULAR_TRANSCRIPTOME } from '../../modules/nf-core/bedtools/getfasta/main'
include { GAWK as COUNTS_TO_BED                       } from '../../modules/nf-core/gawk/main'
include { GNU_SORT as UNIQUE_REGIONS                  } from '../../modules/nf-core/gnu/sort/main'
include { CAT_CAT as CAT_TRANSCRIPTOME                } from '../../modules/nf-core/cat/cat/main'
include { GAWK as MARK_CIRCULAR                       } from '../../modules/nf-core/gawk/main'
include { PSIRC_INDEX                                 } from '../../modules/local/psirc/index/main'
include { PSIRC_QUANT                                 } from '../../modules/local/psirc/quant/main'
include { CUSTOM_TX2GENE                              } from '../../modules/nf-core/custom/tx2gene/main'
include { TXIMETA_TXIMPORT                            } from '../../modules/nf-core/tximeta/tximport/main'
include { COMBINE_QUANTIFICATION as COMBINE_COUNTS    } from '../../modules/local/quantification/combine_quantification/main'
include { COMBINE_QUANTIFICATION as COMBINE_TPM       } from '../../modules/local/quantification/combine_quantification/main'

workflow QUANTIFICATION {
    take:
    ch_gtf
    ch_fasta
    counts
    reads
    circ_annotation

    main:
    ch_versions = Channel.empty()

    LINEAR_TRANSCRIPTOME(ch_gtf, ch_fasta)
    COUNTS_TO_BED(counts.map{counts -> [[id: "counts"], counts]}, [])
    UNIQUE_REGIONS(COUNTS_TO_BED.out.output)

    CIRCULAR_TRANSCRIPTOME(
        UNIQUE_REGIONS.out.sorted.map{meta, bed -> bed},
        ch_fasta.map{meta, fasta -> fasta}
    )

    CAT_TRANSCRIPTOME(
        LINEAR_TRANSCRIPTOME.out.transcriptome.map{meta, fasta -> fasta}.mix(
            CIRCULAR_TRANSCRIPTOME.out.fasta
        ).collect().map(fastas -> [[id: "transcriptome"], fastas])
    )

    MARK_CIRCULAR(CAT_TRANSCRIPTOME.out.file_out, [])

    ch_versions = ch_versions.mix(
        LINEAR_TRANSCRIPTOME.out.versions,
        COUNTS_TO_BED.out.versions,
        UNIQUE_REGIONS.out.versions,
        CIRCULAR_TRANSCRIPTOME.out.versions,
        CAT_TRANSCRIPTOME.out.versions,
        MARK_CIRCULAR.out.versions
    )

    PSIRC_INDEX(MARK_CIRCULAR.out.output)
    PSIRC_QUANT(reads, PSIRC_INDEX.out.index.collect())

    CUSTOM_TX2GENE(
        ch_gtf,
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

    ch_counts = TXIMETA_TXIMPORT.out.counts_gene
    ch_tpm = TXIMETA_TXIMPORT.out.tpm_gene

    COMBINE_COUNTS(
        ch_counts.map{meta, counts -> counts}.collect().map{[[id: "counts"], it]},
        CUSTOM_TX2GENE.out.tx2gene,
        circ_annotation
    )

    COMBINE_TPM(
        ch_tpm.map{meta, tpm -> tpm}.collect().map{[[id: "tpm"], it]},
        CUSTOM_TX2GENE.out.tx2gene,
        circ_annotation
    )

    ch_versions = ch_versions.mix(
        COMBINE_COUNTS.out.versions,
        COMBINE_TPM.out.versions
    )

    emit:
    counts_linear   = COMBINE_COUNTS.out.linear
    counts_circular = COMBINE_COUNTS.out.circular
    counts_combined = COMBINE_COUNTS.out.combined

    tpm_linear      = COMBINE_TPM.out.linear
    tpm_circular    = COMBINE_TPM.out.circular
    tpm_combined    = COMBINE_TPM.out.combined

    versions        = ch_versions
}
