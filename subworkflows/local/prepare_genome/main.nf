include { CUSTOM_GTFFILTER as GTFFILTER                      } from '../../../modules/nf-core/custom/gtffilter'
include { SEQKIT_SPLIT                                       } from '../../../modules/local/seqkit/split'
include { BOWTIE_BUILD                                       } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE2_BUILD                                      } from '../../../modules/nf-core/bowtie2/build'
include { BWA_INDEX                                          } from '../../../modules/nf-core/bwa/index'
include { HISAT2_EXTRACTSPLICESITES                          } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                                       } from '../../../modules/nf-core/hisat2/build'
include { STAR_GENOMEGENERATE                                } from '../../../modules/nf-core/star/genomegenerate'
include { GAWK as CLEAN_FASTA                                } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                                     } from '../../../modules/nf-core/samtools/faidx'
include { UCSC_GTFTOGENEPRED                                 } from '../../../modules/nf-core/ucsc/gtftogenepred'
include { GAWK as CIRCEXPLORER2_REFERENCE                    } from '../../../modules/nf-core/gawk'
include { PSIRC_TRANSCRIPTOME                                } from '../../../modules/local/psirc/transcriptome'
include { PSIRC_INDEX                                        } from '../../../modules/local/psirc/index'
include { starBasedTools                                     } from '../utils_nfcore_circrna_pipeline'
include { CIRCMINER_INDEX                                    } from '../../../modules/local/circminer/index'
include { AGAT_CONVERTSPGFF2GTF as CIRCMINER_STANDARDIZE_GTF } from '../../../modules/nf-core/agat/convertspgff2gtf'
include { CIRCMINER_ORDERGTF                                 } from '../../../modules/local/circminer/ordergtf'

workflow PREPARE_GENOME {
    take:
    ch_fasta
    ch_gtf
    tools
    quantification_tools
    bowtie_index
    bowtie2_index
    bwa_index
    hisat2_index
    star_index
    fli_tools
    circminer_index

    main:
    fli_tools = fli_tools ?: []

    // MapSplice cannot deal with extra field in the fasta headers
    // this removes all additional fields in the headers of the input fasta file
    if (tools.contains('mapsplice')) {
        CLEAN_FASTA(ch_fasta.map { meta, fasta -> [meta, fasta, ''] }, [], false)
        ch_fasta = CLEAN_FASTA.out.output
    }

    GTFFILTER(ch_gtf, ch_fasta)
    ch_gtf = GTFFILTER.out.gtf

    // circMiner needs parent records and transcript-direction exon ordering.
    // Keep this annotation separate: AGAT repairs must not alter other callers' inputs.
    ch_circminer_gtf = channel.empty()
    if (tools.contains('circminer')) {
        CIRCMINER_STANDARDIZE_GTF(ch_gtf)
        CIRCMINER_ORDERGTF(CIRCMINER_STANDARDIZE_GTF.out.output_gtf)
        ch_circminer_gtf = CIRCMINER_ORDERGTF.out.gtf
    }

    UCSC_GTFTOGENEPRED(ch_gtf)

    SEQKIT_SPLIT(ch_fasta)

    if (bowtie_index) {
        ch_bowtie_index = channel.value([[id: "bowtie"], bowtie_index])
    }
    else if (tools.contains('mapsplice')) {
        BOWTIE_BUILD(ch_fasta)
        ch_bowtie_index = BOWTIE_BUILD.out.index
    }
    else {
        ch_bowtie_index = channel.empty()
    }

    if (bowtie2_index) {
        ch_bowtie2_index = channel.value([[id: "bowtie2"], bowtie2_index])
    }
    else if (tools.contains('find_circ')) {
        BOWTIE2_BUILD(ch_fasta)
        ch_bowtie2_index = BOWTIE2_BUILD.out.index
    }
    else {
        ch_bowtie2_index = channel.empty()
    }

    if (bwa_index) {
        ch_bwa_index = channel.value([[id: "bwa"], bwa_index])
    }
    else if (tools.intersect(['ciri2', 'ciri3']).size() > 0 || quantification_tools.intersect(['ciriquant', 'ciri3']).size() > 0) {
        BWA_INDEX(ch_fasta)
        ch_bwa_index = BWA_INDEX.out.index
    }
    else {
        ch_bwa_index = channel.empty()
    }

    ch_hisat2_splice_sites = channel.empty()
    if (hisat2_index) {
        ch_hisat2_index = channel.value([[id: "hisat2"], hisat2_index])
    }
    else if (tools.contains('ciri2') || quantification_tools.contains('ciriquant')) {
        HISAT2_EXTRACTSPLICESITES(ch_gtf)
        ch_hisat2_splice_sites = HISAT2_EXTRACTSPLICESITES.out.txt

        ch_hisat2_reference = ch_fasta
            .combine(ch_gtf)
            .combine(HISAT2_EXTRACTSPLICESITES.out.txt)
            .map { meta_fasta, fasta, _meta_gtf, gtf, _meta_splices, splices ->
                [meta_fasta, fasta, gtf, splices]
            }
        HISAT2_BUILD(ch_hisat2_reference, params.hisat2_build_memory)
        ch_hisat2_index = HISAT2_BUILD.out.index
    }
    else {
        ch_hisat2_index = channel.empty()
    }

    if (star_index) {
        ch_star_index = channel.value([[id: "star"], star_index])
    }
    else if (starBasedTools(tools).size() > 0 || fli_tools.contains('circtools')) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_star_index = STAR_GENOMEGENERATE.out.index
    }
    else {
        ch_star_index = channel.empty()
    }

    SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [meta, fasta, []] }, false)

    // Circexplorer2 reference is needed for annotation
    CIRCEXPLORER2_REFERENCE(UCSC_GTFTOGENEPRED.out.genepred.map { meta, genepred -> [meta, genepred, 'circexplorer2.txt'] }, [], false)
    ch_circexplorer2_reference = CIRCEXPLORER2_REFERENCE.out.output.map { _meta, file -> file }.collect()

    ch_psirc_index = channel.empty()
    if (tools.contains('psirc')) {
        PSIRC_TRANSCRIPTOME(ch_gtf, ch_fasta)

        PSIRC_INDEX(PSIRC_TRANSCRIPTOME.out.transcriptome)
        ch_psirc_index = PSIRC_TRANSCRIPTOME.out.transcriptome
            .join(
                PSIRC_INDEX.out.index.map { meta, a, b -> [meta, [a, b]] }
            )
            .collect()
    }

    if (circminer_index) {
        ch_circminer_index = channel.value([[id: "circminer"], circminer_index])
    }
    else if (tools.contains('circminer')) {
        CIRCMINER_INDEX(ch_fasta)
        ch_circminer_index = CIRCMINER_INDEX.out.index
    }
    else {
        ch_circminer_index = channel.empty()
    }

    emit:
    gtf             = ch_gtf
    faidx           = SAMTOOLS_FAIDX.out.fai
    bowtie_index    = ch_bowtie_index
    bowtie2_index   = ch_bowtie2_index
    bwa_index       = ch_bwa_index
    hisat2_index    = ch_hisat2_index
    star_index      = ch_star_index
    circexplorer2   = ch_circexplorer2_reference
    psirc_index     = ch_psirc_index
    circminer_index = ch_circminer_index
    chromosomes     = SEQKIT_SPLIT.out.split
    splice_sites    = ch_hisat2_splice_sites.collect()
    circminer_gtf   = ch_circminer_gtf
}
