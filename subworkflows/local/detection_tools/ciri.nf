include { BWA_MEM                       } from '../../../modules/nf-core/bwa/mem'
include { CIRIQUANT                     } from '../../../modules/local/ciriquant/ciriquant'
include { CIRI_CIRI2 as CIRI2           } from '../../../modules/local/ciri/ciri2'
include { GAWK as UNIFY                 } from '../../../modules/nf-core/gawk'
include { SEQKIT_FX2TAB                 } from '../../../modules/nf-core/seqkit/fx2tab'
include { CIRI_READLENGTH as READLENGTH } from '../../../modules/local/ciri/readlength'
include { FASTP                         } from '../../../modules/nf-core/fastp'

workflow CIRI {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    ch_bwa_index

    main:
    ch_versions = Channel.empty()

    def cirifull_enabled = params.fli_tools.split(',').collect { it.trim() }.contains('cirifull')

    if (cirifull_enabled) {
        // CIRI-full requires all reads to have the same length

        ch_read1 = ch_reads.map { meta, reads -> [[id: meta.id + '_r1', old_meta: meta, r: 1], reads[0]] }
        ch_read2 = ch_reads.map { meta, reads -> [[id: meta.id + '_r2', old_meta: meta, r: 2], reads[1]] }

        // Get the read lengths
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

        // Determine the 5th percentile of the read lengths
        READLENGTH(ch_reads_len)
        ch_versions = ch_versions.mix(READLENGTH.out.versions)

        ch_fastp = ch_reads
            .join(READLENGTH.out.length)
            .map { meta, reads, length -> [meta + [target_length: length.text.toInteger()], reads] }

        // Trim the reads to the 5th percentile length
        FASTP(ch_fastp, [], false, false, false)
        ch_versions = ch_versions.mix(FASTP.out.versions)

        // These are the new reads to use for the rest of the CIRI pipeline
        ch_reads = FASTP.out.reads
    }

    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    CIRI2(BWA_MEM.out.sam, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRI2.out.versions)

    UNIFY(
        CIRI2.out.txt.map { meta, txt ->
            [meta + [tool: "ciri"], txt]
        },
        [],
        false,
    )
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed                = UNIFY.out.output
    ciri_txt           = CIRI2.out.txt
    ciri_sam           = BWA_MEM.out.sam
    reads_fixed_length = cirifull_enabled ? FASTP.out.reads : Channel.empty()
    versions           = ch_versions
}
