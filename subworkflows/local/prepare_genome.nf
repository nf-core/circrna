include { SEQKIT_SPLIT        } from '../../modules/local/seqkit/split'
include { BOWTIE_BUILD        } from '../../modules/nf-core/bowtie/build'
include { BOWTIE2_BUILD       } from '../../modules/nf-core/bowtie2/build'
include { BWA_INDEX           } from '../../modules/nf-core/bwa/index'
include { HISAT2_EXTRACTSPLICESITES } from '../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD        } from '../../modules/nf-core/hisat2/build'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate'
include { SEGEMEHL_INDEX      } from '../../modules/nf-core/segemehl/index'
include { GAWK as CLEAN_FASTA } from '../../modules/nf-core/gawk'
include { SAMTOOLS_FAIDX      } from '../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {

    take:
    ch_fasta
    ch_gtf

    main:
    ch_versions = Channel.empty()

    // MapSplice cannot deal with extra field in the fasta headers
    // this removes all additional fields in the headers of the input fasta file
    if( params.module.contains('circrna_discovery') && params.tool.contains('mapsplice') ) {
        CLEAN_FASTA(ch_fasta, [])
        ch_fasta = CLEAN_FASTA.out.output

        ch_versions = ch_versions.mix(CLEAN_FASTA.out.versions)
    }

    SEQKIT_SPLIT(ch_fasta)

    BOWTIE_BUILD(ch_fasta.map{ meta, fasta -> fasta })

    BOWTIE2_BUILD(ch_fasta)

    BWA_INDEX (ch_fasta)

    HISAT2_EXTRACTSPLICESITES(ch_gtf)

    HISAT2_BUILD(ch_fasta, ch_gtf, HISAT2_EXTRACTSPLICESITES.out.txt)

    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)

    SEGEMEHL_INDEX(ch_fasta.map{ meta, fasta -> fasta}) // TODO: Add support for meta

    SAMTOOLS_FAIDX(ch_fasta, [[], []])

    // Collect versions
    ch_versions = ch_versions.mix(SEQKIT_SPLIT.out.versions,
                                    BOWTIE_BUILD.out.versions,
                                    BOWTIE2_BUILD.out.versions,
                                    BWA_INDEX.out.versions,
                                    HISAT2_EXTRACTSPLICESITES.out.versions,
                                    HISAT2_BUILD.out.versions,
                                    SEGEMEHL_INDEX.out.versions,
                                    STAR_GENOMEGENERATE.out.versions)

    emit:
    faidx        = SAMTOOLS_FAIDX.out.fai
    bowtie       = params.bowtie   ?: BOWTIE_BUILD.out.index
    segemehl     = params.segemehl ?: SEGEMEHL_INDEX.out.index
    bowtie2      = params.bowtie2  ? Channel.value([[id: "bowtie2"], file(params.bowtie2, checkIfExists: true)]) : BOWTIE2_BUILD.out.index.collect()
    bwa          = params.bwa      ? Channel.value([[id: "bwa"], file(params.bwa, checkIfExists: true)])         : BWA_INDEX.out.index.collect()
    hisat2       = params.hisat2   ? Channel.value([[id: "hisat2"], file(params.hisat2, checkIfExists: true)])   : HISAT2_BUILD.out.index.collect()
    star         = params.star     ? Channel.value([[id: "star"], file(params.star, checkIfExists: true)])       : STAR_GENOMEGENERATE.out.index.collect()
    chromosomes  = SEQKIT_SPLIT.out.split
    splice_sites = HISAT2_EXTRACTSPLICESITES.out.txt.collect()

    versions     = ch_versions
}
