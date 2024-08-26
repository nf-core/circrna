include { CIRI_CIRIFULL as MAIN } from '../../../modules/local/ciri/cirifull'

workflow CIRI_FULL {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    
    main:
    ch_versions = Channel.empty()

    MAIN(ch_reads, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(MAIN.out.versions)
    
    emit:
    versions = ch_versions
}
