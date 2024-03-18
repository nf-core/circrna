include { SEQKIT_SPLIT        } from '../../modules/local/seqkit/split/main'
include { BOWTIE_BUILD        } from '../../modules/nf-core/bowtie/build/main'
include { BOWTIE2_BUILD       } from '../../modules/nf-core/bowtie2/build/main'
include { BWA_INDEX           } from '../../modules/nf-core/bwa/index/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD        } from '../../modules/nf-core/hisat2/build/main'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main'
include { SEGEMEHL_INDEX      } from '../../modules/nf-core/segemehl/index/main'
include { GAWK as CLEAN_FASTA } from '../../modules/nf-core/gawk/main'

workflow PREPARE_GENOME {

    take:
    ch_fasta
    ch_gtf

    main:
    ch_versions = Channel.empty()

    // MapSplice cannot deal with extra field in the fasta headers
    // this removes all additional fields in the headers of the input fasta file
    if( params.tool.contains('mapsplice') && params.module.contains('circrna_discovery') ) {
        CLEAN_FASTA(ch_fasta, [])
        ch_fasta = CLEAN_FASTA.out.output
    }

    SEQKIT_SPLIT(ch_fasta)

    BOWTIE_BUILD(ch_fasta.map{ meta, fasta -> fasta })

    BOWTIE2_BUILD(ch_fasta)

    BWA_INDEX (ch_fasta)

    HISAT2_EXTRACTSPLICESITES(ch_gtf)

    HISAT2_BUILD(ch_fasta, ch_gtf, HISAT2_EXTRACTSPLICESITES.out.txt)

    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)

    SEGEMEHL_INDEX(ch_fasta.map{ meta, fasta -> fasta}) // TODO: Add support for meta

    // Collect versions
    ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
    ch_versions = ch_versions.mix(SEGEMEHL_INDEX.out.versions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    bowtie       = BOWTIE_BUILD.out.index
    bowtie2      = BOWTIE2_BUILD.out.index
    bwa          = BWA_INDEX.out.index
    chromosomes  = SEQKIT_SPLIT.out.split
    hisat2       = HISAT2_BUILD.out.index.collect()
    star         = STAR_GENOMEGENERATE.out.index.collect()
    segemehl     = SEGEMEHL_INDEX.out.index
    splice_sites = HISAT2_EXTRACTSPLICESITES.out.txt.collect()
    versions     = ch_versions
}
