include { BWA_MEM     } from '../../../modules/nf-core/bwa/mem'
include { CIRI_CIRI2  } from '../../../modules/local/ciri/ciri2'
include { CIRI_CIRIAS } from '../../../modules/local/ciri/cirias'

workflow CIRI_FULL {
    take:
    ch_bwa_index
    ch_fasta
    ch_gtf
    ch_reads
    
    main:
    ch_versions = Channel.empty()

    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, false)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    CIRI_CIRI2(BWA_MEM.out.bam, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRI_CIRI2.out.versions)

    CIRI_CIRIAS(BWA_MEM.out.bam.join(CIRI_CIRI2.out.ciri), ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRI_CIRIAS.out.versions)
    
    emit:
    versions = ch_versions
}