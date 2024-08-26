include { CIRI_CIRIFULL as MAIN } from '../../../modules/local/ciri/cirifull'

workflow CIRI_FULL {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    
    main:
    ch_versions = Channel.empty()

    // TODO: Make sure all reads have the same length before running CIRI-full
    // https://www.biostars.org/p/332512/

    MAIN(ch_reads, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(MAIN.out.versions)
    
    emit:
    versions = ch_versions
}
