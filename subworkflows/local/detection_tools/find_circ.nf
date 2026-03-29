include { BOWTIE2_ALIGN as ALIGN       } from '../../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_VIEW                } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX               } from '../../../modules/nf-core/samtools/index'
include { FIND_CIRC_ANCHORS as ANCHORS } from '../../../modules/local/find_circ/anchors'
include { FIND_CIRC as MAIN            } from '../../../modules/local/find_circ/find_circ'
include { BIOAWK as EXTRACT_READS      } from '../../../modules/nf-core/bioawk'
include { CSVTK_SUMMARY as GROUP_READS } from '../../../modules/local/csvtk/summary'
include { CSVTK_JOIN as JOIN_READS     } from '../../../modules/nf-core/csvtk/join'
include { GAWK as UNIFY                } from '../../../modules/nf-core/gawk'

workflow FIND_CIRC {
    take:
    reads
    bowtie2_index
    ch_fasta

    main:
    ch_versions = Channel.empty()

    ALIGN( reads, bowtie2_index, ch_fasta, false, true )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    SAMTOOLS_INDEX( ALIGN.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    SAMTOOLS_VIEW( ALIGN.out.bam.join( SAMTOOLS_INDEX.out.bai ), ch_fasta, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    ANCHORS( SAMTOOLS_VIEW.out.bam )
    ch_versions = ch_versions.mix(ANCHORS.out.versions)

    MAIN( ANCHORS.out.anchors, bowtie2_index, ch_fasta.map{ _meta, fasta -> fasta } )
    ch_versions = ch_versions.mix(MAIN.out.versions)

    EXTRACT_READS( MAIN.out.reads )
    ch_versions = ch_versions.mix(EXTRACT_READS.out.versions)

    GROUP_READS( EXTRACT_READS.out.output )
    ch_versions = ch_versions.mix(GROUP_READS.out.versions)

    JOIN_READS( MAIN.out.bed.join(GROUP_READS.out.csv)
        .map{ meta, bed, _reads -> [ meta, [bed, _reads]] }
    )
    ch_versions = ch_versions.mix(JOIN_READS.out.versions)

    UNIFY( JOIN_READS.out.csv.map{ meta, bed ->
        [ meta + [tool: "find_circ"], bed ] }, [], false )
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}
