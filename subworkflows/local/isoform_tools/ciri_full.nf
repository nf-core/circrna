include { BWA_MEM     } from '../../../modules/nf-core/bwa/mem'
include { CIRI_CIRI2  } from '../../../modules/local/ciri/ciri2'
include { CIRI_CIRIAS } from '../../../modules/local/ciri/cirias'

workflow CIRI_FULL {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    ch_bwa_index
    ch_bsjs

    main:
    ch_versions = Channel.empty()

    // TODO: Make sure all reads have the same length before running BWA-MEM
    // https://www.biostars.org/p/332512/
    // Could utilize BIOAWK for this
    // https://medium.com/@odettegeorge/bioawk-33a956e483ec

    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, false)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    // TODO: Replcase CIRI2 with a process that transforms the BED file in ch_bsjs
    // into a compatible format. This ensures that we can control what BSJs should
    // be used as base for FLI detection.

    CIRI_CIRI2(BWA_MEM.out.bam, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRI_CIRI2.out.versions)

    CIRI_CIRIAS(BWA_MEM.out.bam.join(CIRI_CIRI2.out.ciri), ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRI_CIRIAS.out.versions)

    // TODO: Implement remaining steps described in
    // https://ciri-cookbook.readthedocs.io/en/latest/CIRI-full.html#running-ciri-full-step-by-step
    
    emit:
    versions = ch_versions
}
