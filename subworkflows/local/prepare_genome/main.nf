include { CUSTOM_GTFFILTER          as GTFFILTER                 } from '../../../modules/nf-core/custom/gtffilter'
include { SEQKIT_SPLIT                                           } from '../../../modules/local/seqkit/split'
include { BOWTIE_BUILD                                           } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE2_BUILD                                          } from '../../../modules/nf-core/bowtie2/build'
include { BWA_INDEX                                              } from '../../../modules/nf-core/bwa/index'
include { HISAT2_EXTRACTSPLICESITES                              } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                                           } from '../../../modules/nf-core/hisat2/build'
include { STAR_GENOMEGENERATE                                    } from '../../../modules/nf-core/star/genomegenerate'
include { GAWK                      as CLEAN_FASTA               } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                                         } from '../../../modules/nf-core/samtools/faidx'
include { UCSC_GTFTOGENEPRED                                     } from '../../../modules/nf-core/ucsc/gtftogenepred'
include { GAWK                      as CIRCEXPLORER2_REFERENCE   } from '../../../modules/nf-core/gawk'
include { GFFREAD                                                } from '../../../modules/nf-core/gffread'
include { PSIRC_TRANSCRIPTOME                                    } from '../../../modules/local/psirc/transcriptome'
include { PSIRC_INDEX                                            } from '../../../modules/local/psirc/index'
include { starBasedTools                                         } from '../utils_nfcore_circrna_pipeline'
include { CIRCMINER_INDEX                                        } from '../../../modules/local/circminer/index'
include { AGAT_CONVERTSPGFF2GTF     as CIRCMINER_STANDARDIZE_GTF } from '../../../modules/nf-core/agat/convertspgff2gtf'
include { CIRCMINER_ORDERGTF                                     } from '../../../modules/local/circminer/ordergtf'

workflow PREPARE_GENOME {
    take:
    ch_fasta
    ch_gtf
    tools
    quantification_tools
    bowtie
    bowtie2
    bwa
    hisat2
    star
    fli_tools
    circminer

    main:
    fli_tools = fli_tools ?: []
    ch_versions = channel.empty()

    // MapSplice cannot deal with extra field in the fasta headers
    // this removes all additional fields in the headers of the input fasta file
    if (tools.contains('mapsplice')) {
        CLEAN_FASTA(ch_fasta.map { meta, fasta -> [meta, fasta, ''] }, [], false)
        ch_fasta = CLEAN_FASTA.out.output
    }

    GTFFILTER(ch_gtf, ch_fasta)
    ch_gtf = GTFFILTER.out.gtf
    ch_versions = ch_versions.mix(GTFFILTER.out.versions)

    // circMiner needs parent records and transcript-direction exon ordering.
    // Keep this annotation separate: AGAT repairs must not alter other callers' inputs.
    ch_circminer_gtf = channel.empty()
    if (tools.contains('circminer')) {
        CIRCMINER_STANDARDIZE_GTF(ch_gtf)
        CIRCMINER_ORDERGTF(CIRCMINER_STANDARDIZE_GTF.out.output_gtf)
        ch_circminer_gtf = CIRCMINER_ORDERGTF.out.gtf
        ch_versions = ch_versions.mix(CIRCMINER_ORDERGTF.out.versions)
    }

    UCSC_GTFTOGENEPRED(ch_gtf)

    SEQKIT_SPLIT(ch_fasta)
    ch_versions = ch_versions.mix(SEQKIT_SPLIT.out.versions)

    if (bowtie) {
        ch_bowtie = channel.value([[id: "bowtie"], bowtie])
    }
    else if (tools.contains('mapsplice')) {
        BOWTIE_BUILD(ch_fasta)
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
        ch_bowtie = BOWTIE_BUILD.out.index
    }
    else {
        ch_bowtie = channel.empty()
    }

    if (bowtie2) {
        ch_bowtie2 = channel.value([[id: "bowtie2"], bowtie2])
    }
    else if (tools.contains('find_circ')) {
        BOWTIE2_BUILD(ch_fasta)
        ch_bowtie2 = BOWTIE2_BUILD.out.index
    }
    else {
        ch_bowtie2 = channel.empty()
    }

    if (bwa) {
        ch_bwa = channel.value([[id: "bwa"], bwa])
    }
    else if (tools.contains('ciri') || tools.contains('ciri3') || quantification_tools.contains('ciriquant')) {
        BWA_INDEX(ch_fasta)
        ch_bwa = BWA_INDEX.out.index
    }
    else {
        ch_bwa = channel.empty()
    }

    ch_hisat2_splice_sites = channel.empty()
    if (hisat2) {
        ch_hisat2 = channel.value([[id: "hisat2"], hisat2])
    }
    else if (tools.contains('ciri') || quantification_tools.contains('ciriquant')) {
        HISAT2_EXTRACTSPLICESITES(ch_gtf)
        ch_hisat2_splice_sites = HISAT2_EXTRACTSPLICESITES.out.txt

        HISAT2_BUILD(ch_fasta, ch_gtf, HISAT2_EXTRACTSPLICESITES.out.txt)
        ch_hisat2 = HISAT2_BUILD.out.index
    }
    else {
        ch_hisat2 = channel.empty()
    }

    if (star) {
        ch_star = channel.value([[id: "star"], star])
    }
    else if (starBasedTools(tools).size() > 0 || fli_tools.contains('circtools')) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_star = STAR_GENOMEGENERATE.out.index
    }
    else {
        ch_star = channel.empty()
    }

    SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [meta, fasta, []] }, false)

    // Circexplorer2 reference is needed for annotation
    CIRCEXPLORER2_REFERENCE(UCSC_GTFTOGENEPRED.out.genepred.map { meta, genepred -> [meta, genepred, 'circexplorer2.txt'] }, [], false)
    ch_circexplorer2_reference = CIRCEXPLORER2_REFERENCE.out.output.map { _meta, file -> file }.collect()

    ch_psirc_index = channel.empty()
    if (tools.contains('psirc')) {
        PSIRC_TRANSCRIPTOME(ch_gtf, ch_fasta)
        ch_versions = ch_versions.mix(PSIRC_TRANSCRIPTOME.out.versions)

        PSIRC_INDEX(PSIRC_TRANSCRIPTOME.out.transcriptome)
        ch_versions = ch_versions.mix(PSIRC_INDEX.out.versions)
        ch_psirc_index = PSIRC_TRANSCRIPTOME.out.transcriptome.join(
            PSIRC_INDEX.out.index.map { meta, a, b -> [meta, [a, b]]}
        ).collect()
    }

    if (circminer) {
        ch_circminer = channel.value([[id: "circminer"], circminer])
    } else if (tools.contains('circminer')) {
        CIRCMINER_INDEX(ch_fasta)
        ch_versions = ch_versions.mix(CIRCMINER_INDEX.out.versions)
        ch_circminer = CIRCMINER_INDEX.out.index
    } else {
        ch_circminer = channel.empty()
    }

    emit:
    gtf           = ch_gtf
    faidx         = SAMTOOLS_FAIDX.out.fai
    bowtie        = ch_bowtie
    bowtie2       = ch_bowtie2
    bwa           = ch_bwa
    hisat2        = ch_hisat2
    star          = ch_star
    circexplorer2 = ch_circexplorer2_reference
    psirc         = ch_psirc_index
    circminer     = ch_circminer
    chromosomes   = SEQKIT_SPLIT.out.split
    splice_sites  = ch_hisat2_splice_sites.collect()
    versions      = ch_versions
    circminer_gtf = ch_circminer_gtf
}
