include { CIRIQUANT as MAIN                } from '../../../modules/local/ciriquant'
include { PYGTFTK_TABULATE as EXTRACT_CIRC } from '../../../modules/local/pygtftk/tabulate'
include { GAWK as EXTRACT_GENES            } from '../../../modules/nf-core/gawk'

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

    EXTRACT_GENES( MAIN.out.gene, [] )
    ch_versions = ch_versions.mix(EXTRACT_GENES.out.versions)
    
    emit:


    versions = ch_versions
}
