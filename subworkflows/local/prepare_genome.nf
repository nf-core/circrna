include { CUSTOM_GTFFILTER as GTFFILTER   } from '../../modules/nf-core/custom/gtffilter'
include { SEQKIT_SPLIT                    } from '../../modules/local/seqkit/split'
include { BOWTIE_BUILD                    } from '../../modules/nf-core/bowtie/build'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index'
include { HISAT2_EXTRACTSPLICESITES       } from '../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                    } from '../../modules/nf-core/hisat2/build'
include { STAR_GENOMEGENERATE             } from '../../modules/nf-core/star/genomegenerate'
include { GAWK as CLEAN_FASTA             } from '../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx'
include { UCSC_GTFTOGENEPRED              } from '../../modules/nf-core/ucsc/gtftogenepred'
include { GAWK as CIRCEXPLORER2_REFERENCE } from '../../modules/nf-core/gawk'
include { GFFREAD                         } from '../../modules/nf-core/gffread'
include { PSIRC_TRANSCRIPTOME             } from '../../modules/local/psirc/transcriptome'
include { PSIRC_INDEX                     } from '../../modules/local/psirc/index'

workflow PREPARE_GENOME {
    take:
    ch_fasta
    ch_gtf

    main:
    ch_versions = Channel.empty()

    detection_tools = params.tools.split(',').collect { it.trim().toLowerCase() }

    // MapSplice cannot deal with extra field in the fasta headers
    // this removes all additional fields in the headers of the input fasta file
    if (detection_tools.contains('mapsplice')) {
        CLEAN_FASTA(ch_fasta, [], false)
        ch_fasta = CLEAN_FASTA.out.output

        ch_versions = ch_versions.mix(CLEAN_FASTA.out.versions)
    }

    GTFFILTER(ch_gtf, ch_fasta)
    ch_gtf = GTFFILTER.out.gtf
    ch_versions = ch_versions.mix(GTFFILTER.out.versions)

    UCSC_GTFTOGENEPRED(ch_gtf)
    ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)

    SEQKIT_SPLIT(ch_fasta)
    ch_versions = ch_versions.mix(SEQKIT_SPLIT.out.versions)

    if (params.bowtie) {
        ch_bowtie = Channel.value([[id: "bowtie"], file(params.bowtie, checkIfExists: true)])
    }
    else if (detection_tools.contains('mapsplice')) {
        BOWTIE_BUILD(ch_fasta)
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
        ch_bowtie = BOWTIE_BUILD.out.index
    }
    else {
        ch_bowtie = Channel.empty()
    }

    if (params.bowtie2) {
        ch_bowtie2 = Channel.value([[id: "bowtie2"], file(params.bowtie2, checkIfExists: true)])
    }
    else if (detection_tools.contains('find_circ')) {
        BOWTIE2_BUILD(ch_fasta)
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        ch_bowtie2 = BOWTIE2_BUILD.out.index
    }
    else {
        ch_bowtie2 = Channel.empty()
    }

    if (params.bwa) {
        ch_bwa = Channel.value([[id: "bwa"], file(params.bwa, checkIfExists: true)])
    }
    else if (detection_tools.contains('ciri')) {
        BWA_INDEX(ch_fasta)
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        ch_bwa = BWA_INDEX.out.index
    }
    else {
        ch_bwa = Channel.empty()
    }

    ch_hisat2_splice_sites = Channel.empty()
    if (params.hisat2) {
        ch_hisat2 = Channel.value([[id: "hisat2"], file(params.hisat2, checkIfExists: true)])
    }
    else if (detection_tools.contains('ciri')) {
        HISAT2_EXTRACTSPLICESITES(ch_gtf)
        ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        ch_hisat2_splice_sites = HISAT2_EXTRACTSPLICESITES.out.txt

        HISAT2_BUILD(ch_fasta, ch_gtf, HISAT2_EXTRACTSPLICESITES.out.txt)
        ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
        ch_hisat2 = HISAT2_BUILD.out.index
    }
    else {
        ch_hisat2 = Channel.empty()
    }

    if (params.star) {
        ch_star = Channel.value([[id: "star"], file(params.star, checkIfExists: true)])
    }
    else if (detection_tools.intersect(['circexplorer2', 'circrna_finder', 'dcc', 'mapsplice']).size() > 0) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        ch_star = STAR_GENOMEGENERATE.out.index
    }
    else {
        ch_star = Channel.empty()
    }

    SAMTOOLS_FAIDX(ch_fasta, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Circexplorer2 reference is needed for annotation
    CIRCEXPLORER2_REFERENCE(UCSC_GTFTOGENEPRED.out.genepred, [], false)
    ch_circexplorer2_reference = CIRCEXPLORER2_REFERENCE.out.output.map { _meta, file -> file }.collect()
    ch_versions = ch_versions.mix(CIRCEXPLORER2_REFERENCE.out.versions)

    ch_psirc_index = Channel.empty()
    if (detection_tools.contains('psirc')) {
        PSIRC_TRANSCRIPTOME(ch_gtf, ch_fasta)
        ch_versions = ch_versions.mix(PSIRC_TRANSCRIPTOME.out.versions)

        PSIRC_INDEX(PSIRC_TRANSCRIPTOME.out.transcriptome)
        ch_versions = ch_versions.mix(PSIRC_INDEX.out.versions)
        ch_psirc_index = PSIRC_TRANSCRIPTOME.out.transcriptome.join(
            PSIRC_INDEX.out.index.map { meta, a, b -> [meta, [a, b]]}
        ).collect()
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
    chromosomes   = SEQKIT_SPLIT.out.split
    splice_sites  = ch_hisat2_splice_sites.collect()
    versions      = ch_versions
}
