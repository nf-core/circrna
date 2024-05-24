include { BOWTIE2_ALIGN as ALIGN       } from '../../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_VIEW                } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX               } from '../../../modules/nf-core/samtools/index'
include { FIND_CIRC_ANCHORS as ANCHORS } from '../../../modules/local/find_circ/anchors'
include { FIND_CIRC as RUN             } from '../../../modules/local/find_circ/find_circ'
include { GAWK as UNIFY                } from '../../../modules/nf-core/gawk'

workflow FIND_CIRC {
    take:
    reads
    bowtie2_index
    ch_fasta
    
    main:
    ch_versions = Channel.empty()

    ALIGN( reads, bowtie2_index, ch_fasta, false, true )
    SAMTOOLS_INDEX( ALIGN.out.bam )
    SAMTOOLS_VIEW( ALIGN.out.bam.join( SAMTOOLS_INDEX.out.bai ), ch_fasta, [] )
    ANCHORS( SAMTOOLS_VIEW.out.bam )
    RUN( ANCHORS.out.anchors, bowtie2_index, ch_fasta.map{ meta, fasta -> fasta } )
    UNIFY( RUN.out.bed.map{ meta, bed ->
        [ meta + [tool: "find_circ"], bed ] }, [] )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_versions = ch_versions.mix(ANCHORS.out.versions)
    ch_versions = ch_versions.mix(RUN.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}