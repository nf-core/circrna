include { BWA_MEM               } from '../../../modules/nf-core/bwa/mem'
include { CIRI_CIRI2 as CIRI2   } from '../../../modules/local/ciri/ciri2'
include { CIRI_CIRIAS as CIRIAS } from '../../../modules/local/ciri/cirias'

include { FILTLONG              } from '../../../modules/nf-core/filtlong'
include { SEQKIT_STATS          } from '../../../modules/nf-core/seqkit/stats'
include { FASTP                 } from '../../../modules/nf-core/fastp'
include { CIRIFULL_RO1          } from '../../../modules/local/cirifull/ro1'

workflow CIRI {
    take:
    ch_reads
    ch_fasta
    ch_gtf
    ch_bwa_index

    main:
    ch_versions = Channel.empty()

    // BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, true)
    // ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    // CIRI2(BWA_MEM.out.bam, ch_fasta, ch_gtf)
    // CIRIAS(BWA_MEM.out.bam, ch_fasta, ch_gtf)

    ch_read1 = ch_reads.map{ meta, reads -> [[id: meta.id + ':r1', old_meta: meta, r: 1], reads[0]] }
    ch_read2 = ch_reads.map{ meta, reads -> [[id: meta.id + ':r2', old_meta: meta, r: 2], reads[1]] }

    FILTLONG(ch_read1.mix(ch_read2).map{ meta, reads -> [meta, [], reads]})
    ch_versions = ch_versions.mix(FILTLONG.out.versions)

    SEQKIT_STATS(FILTLONG.out.reads)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)

    ch_min_len = SEQKIT_STATS.out.stats.map{ meta, stats -> {
        def lines = stats.readLines()
        def header = lines[0].split('\t')
        def columnIndex = header.findIndexOf { it == 'min_len' }
        if (columnIndex == -1) {
            error "Column 'min_len' not found in file $stats"
        }
        def minValue = lines.drop(1).collect { line ->
            line.split('\t')[columnIndex].toInteger()
        }.min()
        [meta, minValue]
    }}

    ch_reads_fastp = FILTLONG.out.reads.join(ch_min_len)
        .map{ meta, reads, min_len -> [meta + [min_len: min_len, single_end: true], reads]}

    FASTP(ch_reads_fastp, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    ch_read1_filt = FASTP.out.reads
        .filter{ meta, _reads -> meta.r == 1 }
        .map{ meta, reads -> [meta.old_meta, reads] }
    ch_read2_filt = FASTP.out.reads
        .filter{ meta, _reads -> meta.r == 2 }
        .map{ meta, reads -> [meta.old_meta, reads] }

    ch_reads_filt = ch_read1_filt
        .join(ch_read2_filt)
        .map{ meta, r1, r2 -> [meta, [r1, r2]] }

    CIRIFULL_RO1(ch_reads_filt)
    //ch_versions = ch_versions.mix(CIRIFULL_RO1.out.versions)


    emit:
    versions = ch_versions
}