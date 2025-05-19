include { CIRI_CIRIAS as CIRIAS    } from '../../../modules/local/ciri/cirias'
include { CIRIFULL_RO1 as RO1      } from '../../../modules/local/cirifull/ro1'
include { BWA_MEM                  } from '../../../modules/nf-core/bwa/mem'
include { CIRIFULL_RO2 as RO2      } from '../../../modules/local/cirifull/ro2'
include { CIRIFULL_MERGE as MERGE  } from '../../../modules/local/cirifull/merge'
include { CIRI_CIRIVIS as CIRI_VIS } from '../../../modules/local/ciri/cirivis'

workflow CIRIFULL {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    ch_bwa_index
    ch_ciri_txt
    ch_ciri_sam

    main:
    ch_versions = Channel.empty()
    CIRIAS(ch_ciri_txt.join(ch_ciri_sam), ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRIAS.out.versions)

    RO1(ch_reads)
    ch_versions = ch_versions.mix(RO1.out.versions)

    BWA_MEM(RO1.out.fastq, ch_bwa_index, ch_fasta, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    RO2(BWA_MEM.out.sam.map { meta, bam -> [meta, bam, meta.target_length] }, ch_fasta)
    ch_versions = ch_versions.mix(RO2.out.versions)

    ch_merge = ch_ciri_txt.join(CIRIAS.out.list).join(RO2.out.list)

    MERGE(ch_merge, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(MERGE.out.versions)

    ch_grouped = MERGE.out.anno
        .join(CIRIAS.out.library_length)
        .map { _meta, anno, library_length -> [[id: 'cirifull'], anno, library_length] }
        .groupTuple()
        .map { meta, anno, library_length -> [meta, anno, library_length, []] }

    // CIRI_VIS(ch_grouped, ch_fasta)
    // ch_versions = ch_versions.mix(CIRI_VIS.out.versions)

    emit:
    versions = ch_versions
}
