include { EXTRACT_TRANSCRIPTOME                    } from '../../modules/local/quantification/transcriptome/main'
include { BEDTOOLS_GETFASTA as TRANSCRIPTOME_FASTA } from '../../modules/nf-core/bedtools/getfasta/main'
include { GAWK as COUNTS_TO_BED                    } from '../../modules/nf-core/gawk/main'
include { CAT_CAT as CAT_TRANSCRIPTOME             } from '../../modules/nf-core/cat/cat/main'
include { GAWK as MARK_CIRCULAR                    } from '../../modules/nf-core/gawk/main'
include { PSIRC_INDEX                              } from '../../modules/local/psirc/index/main'
include { PSIRC_QUANT                              } from '../../modules/local/psirc/quant/main'
include { PSIRC_COMBINE                            } from '../../modules/local/psirc/combine/main'

workflow QUANTIFICATION {
    take:
        ch_gtf
        ch_fasta
        counts
        reads

    main:
        EXTRACT_TRANSCRIPTOME(ch_gtf)
        COUNTS_TO_BED(counts.map{counts -> [[id: "counts"], counts]}, [])

        CAT_TRANSCRIPTOME(
            EXTRACT_TRANSCRIPTOME.out.bed.map{meta, bed -> bed}.mix(
                COUNTS_TO_BED.out.output.map{meta, bed -> bed}
            ).collect().map(beds -> [[id: "transcriptome"], beds])
        )

        TRANSCRIPTOME_FASTA(
            CAT_TRANSCRIPTOME.out.file_out.map{meta, bed -> bed},
            ch_fasta.map{meta, fasta -> fasta}
        )
        ch_transcriptome_fasta = TRANSCRIPTOME_FASTA.out.fasta.map{fasta -> [[id: "transcriptome"], fasta]}

        MARK_CIRCULAR(ch_transcriptome_fasta, [])

        PSIRC_INDEX(MARK_CIRCULAR.out.output)
        PSIRC_QUANT(reads, PSIRC_INDEX.out.collect())
        PSIRC_COMBINE( PSIRC_QUANT.out.abundance_tsv.map{ meta, data -> data }.collect() )
}
