include { TRANSCRIPTOME as LINEAR_TRANSCRIPTOME    } from '../../modules/local/quantification/transcriptome/main'
include { BEDTOOLS_GETFASTA as CIRCULAR_TRANSCRIPTOME } from '../../modules/nf-core/bedtools/getfasta/main'
include { GAWK as COUNTS_TO_BED                    } from '../../modules/nf-core/gawk/main'
include { CAT_CAT as CAT_TRANSCRIPTOME             } from '../../modules/nf-core/cat/cat/main'
include { GAWK as MARK_CIRCULAR                    } from '../../modules/nf-core/gawk/main'
include { PSIRC_INDEX                              } from '../../modules/local/psirc/index/main'
include { PSIRC_QUANT                              } from '../../modules/local/psirc/quant/main'
include { PSIRC_COMBINE                            } from '../../modules/local/psirc/combine/main'
include { CUSTOM_TX2GENE                           } from '../../modules/nf-core/custom/tx2gene/main'
include { TXIMETA_TXIMPORT                         } from '../../modules/nf-core/tximeta/tximport/main'

workflow QUANTIFICATION {
    take:
        ch_gtf
        ch_fasta
        counts
        reads

    main:
        LINEAR_TRANSCRIPTOME(ch_gtf, ch_fasta)
        COUNTS_TO_BED(counts.map{counts -> [[id: "counts"], counts]}, [])

        CIRCULAR_TRANSCRIPTOME(
            COUNTS_TO_BED.out.output.map{meta, bed -> bed},
            ch_fasta.map{meta, fasta -> fasta}
        )

        CAT_TRANSCRIPTOME(
            LINEAR_TRANSCRIPTOME.out.transcriptome.map{meta, fasta -> fasta}.mix(
                CIRCULAR_TRANSCRIPTOME.out.fasta
            ).collect().map(fastas -> [[id: "transcriptome"], fastas])
        )

        MARK_CIRCULAR(CAT_TRANSCRIPTOME.out.file_out, [])

        PSIRC_INDEX(MARK_CIRCULAR.out.output)
        PSIRC_QUANT(reads, PSIRC_INDEX.out.collect())
        // PSIRC_COMBINE( PSIRC_QUANT.out.abundance_tsv.map{ meta, data -> data }.collect() )

        CUSTOM_TX2GENE(
            ch_gtf,
            PSIRC_QUANT.out.directory,
            "kallisto",
            "gene_id",
            "gene_name"
        )
}
