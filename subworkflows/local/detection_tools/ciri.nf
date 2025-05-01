include { BWA_MEM as BWA_MEM_1          } from '../../../modules/nf-core/bwa/mem'
include { CIRI_CIRI2 as CIRI2           } from '../../../modules/local/ciri/ciri2'
include { CIRI_CIRIAS as CIRIAS         } from '../../../modules/local/ciri/cirias'

include { SEQKIT_FX2TAB                 } from '../../../modules/nf-core/seqkit/fx2tab'
include { CIRI_READLENGTH as READLENGTH } from '../../../modules/local/ciri/readlength'
include { FASTP                         } from '../../../modules/nf-core/fastp'
include { CIRIFULL_RO1                  } from '../../../modules/local/cirifull/ro1'
include { BWA_MEM as BWA_MEM_2          } from '../../../modules/nf-core/bwa/mem'
include { SAMTOOLS_VIEW as BAM_TO_SAM   } from '../../../modules/nf-core/samtools/view'
include { CIRIFULL_RO2                  } from '../../../modules/local/cirifull/ro2'

workflow CIRI {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    ch_bwa_index

    main:
    ch_versions = Channel.empty()

    BWA_MEM_1(ch_reads, ch_bwa_index, ch_fasta, true)
    ch_versions = ch_versions.mix(BWA_MEM_1.out.versions)

    CIRI2(BWA_MEM_1.out.bam, ch_fasta, ch_gtf)
    CIRIAS(BWA_MEM_1.out.bam, ch_fasta, ch_gtf)

    ch_read1 = ch_reads.map { meta, reads -> [[id: meta.id + '_r1', old_meta: meta, r: 1], reads[0]] }
    ch_read2 = ch_reads.map { meta, reads -> [[id: meta.id + '_r2', old_meta: meta, r: 2], reads[1]] }

    SEQKIT_FX2TAB(ch_read1.mix(ch_read2))
    ch_versions = ch_versions.mix(SEQKIT_FX2TAB.out.versions)

    ch_read1_len = SEQKIT_FX2TAB.out.text
        .filter { meta, _lengths -> meta.r == 1 }
        .map { meta, lengths -> [meta.old_meta, lengths] }
    ch_read2_len = SEQKIT_FX2TAB.out.text
        .filter { meta, _lengths -> meta.r == 2 }
        .map { meta, lengths -> [meta.old_meta, lengths] }

    ch_reads_len = ch_read1_len
        .join(ch_read2_len)
        .map { meta, r1, r2 -> [meta, [r1, r2]] }

    READLENGTH(ch_reads_len)
    ch_versions = ch_versions.mix(READLENGTH.out.versions)

    ch_fastp = ch_reads.join(READLENGTH.out.length)
        .map { meta, reads, length -> [meta + [target_length: length.text.toInteger()], reads] }

    FASTP(ch_fastp, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    CIRIFULL_RO1(FASTP.out.reads)
    ch_versions = ch_versions.mix(CIRIFULL_RO1.out.versions)

    BWA_MEM_2(CIRIFULL_RO1.out.fastq, ch_bwa_index, ch_fasta, true)
    ch_versions = ch_versions.mix(BWA_MEM_2.out.versions)

    BAM_TO_SAM(BWA_MEM_2.out.bam.map{ meta, bam -> [meta, bam, []]}, ch_fasta,  [])
    ch_versions = ch_versions.mix(BAM_TO_SAM.out.versions)

    // RO2 has issues with reading the SAM file

    // CIRIFULL_RO2(BAM_TO_SAM.out.sam.map{ meta, bam -> [meta, bam, meta.target_length]}, ch_fasta)
    // ch_versions = ch_versions.mix(CIRIFULL_RO2.out.versions)

    emit:
    versions = ch_versions
}
