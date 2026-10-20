include { SEQKIT_FX2TAB                 } from '../../../modules/nf-core/seqkit/fx2tab'
include { CIRI_READLENGTH as READLENGTH } from '../../../modules/local/ciri/readlength'
include { FASTP                         } from '../../../modules/nf-core/fastp'

workflow CIRI_READS {
    take:
    ch_reads

    main:

    // CIRI-full requires every read to have the same length.
    ch_read1 = ch_reads.map { meta, reads -> [[id: meta.id + '_r1', old_meta: meta, r: 1], reads[0]] }
    ch_read2 = ch_reads.map { meta, reads -> [[id: meta.id + '_r2', old_meta: meta, r: 2], reads[1]] }

    SEQKIT_FX2TAB(ch_read1.mix(ch_read2), 'txt')

    ch_read1_len = SEQKIT_FX2TAB.out.text
        .filter { meta, _lengths -> meta.r == 1 }
        .map { meta, lengths -> [meta.old_meta, lengths] }
    ch_read2_len = SEQKIT_FX2TAB.out.text
        .filter { meta, _lengths -> meta.r == 2 }
        .map { meta, lengths -> [meta.old_meta, lengths] }

    ch_reads_len = ch_read1_len
        .join(ch_read2_len)
        .map { meta, r1, r2 -> [meta, [r1, r2]] }

    // The fifth percentile preserves the historical CIRI-full trimming behaviour.
    READLENGTH(ch_reads_len)

    ch_fastp = ch_reads
        .join(READLENGTH.out.length)
        .map { meta, reads, length -> [meta + [target_length: length.text.toInteger()], reads] }

    FASTP(ch_fastp.map { meta, reads -> [meta, reads, []] }, false, false, false)

    emit:
    reads    = FASTP.out.reads
}
