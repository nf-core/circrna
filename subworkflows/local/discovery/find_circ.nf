include { BOWTIE2_ALIGN as ALIGN       } from '../../../modules/nf-core/bowtie2/align'
include { SAMTOOLS_VIEW                } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX               } from '../../../modules/nf-core/samtools/index'
include { FIND_CIRC_ANCHORS as ANCHORS } from '../../../modules/local/find_circ/anchors'
include { FIND_CIRC as RUN             } from '../../../modules/local/find_circ/find_circ'
include { FIND_CIRC_FILTER as FILTER   } from '../../../modules/local/find_circ/filter'


workflow FIND_CIRC {
    take:
    reads
    bowtie2_index
    ch_fasta
    bsj_reads
    
    main:
    ch_versions = Channel.empty()

    ALIGN( reads, bowtie2_index, ch_fasta, false, true )
    SAMTOOLS_INDEX( ALIGN.out.bam )
    SAMTOOLS_VIEW( ALIGN.out.bam.join( SAMTOOLS_INDEX.out.bai ), ch_fasta, [] )
    ANCHORS( SAMTOOLS_VIEW.out.bam )
    RUN( ANCHORS.out.anchors, bowtie2_index, ch_fasta.map{ meta, fasta -> fasta } )
    FILTER( RUN.out.bed.map{ meta, bed ->
        [ meta + [tool: "find_circ"], bed ] }, bsj_reads )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_versions = ch_versions.mix(ANCHORS.out.versions)
    ch_versions = ch_versions.mix(RUN.out.versions)
    ch_versions = ch_versions.mix(FILTER.out.versions)

    emit:
    matrix = FILTER.out.matrix
    results = FILTER.out.results

    versions = ch_versions
}