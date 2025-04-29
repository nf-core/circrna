include { BWA_MEM             } from '../../../modules/nf-core/bwa/mem'
include { CIRI_CIRI2 as CIRI2 } from '../../../modules/local/ciri/ciri2'

workflow CIRI {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    ch_bwa_index

    main:
    ch_versions = Channel.empty()

    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    CIRI2(BWA_MEM.out.bam, ch_fasta, ch_gtf)

    emit:
    versions = ch_versions
}