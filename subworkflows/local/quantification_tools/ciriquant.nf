include { CIRIQUANT as MAIN      } from '../../../modules/local/ciriquant'
include { BIOAWK as EXTRACT_CIRC } from '../../../modules/nf-core/bioawk'

workflow CIRIQUANT {
    take:
    reads
    ch_bed
    ch_gtf
    ch_fasta
    bwa_index
    hisat2_index

    main:
    ch_versions = Channel.empty()

    MAIN( reads, ch_bed, ch_gtf, ch_fasta, bwa_index, hisat2_index )
    ch_versions = ch_versions.mix(MAIN.out.versions)

    EXTRACT_CIRC( MAIN.out.gtf )
    ch_versions = ch_versions.mix(EXTRACT_CIRC.out.versions)
    
    emit:


    versions = ch_versions
}