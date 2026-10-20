include { CIRIFULL } from '../fli_tools/cirifull'
include { PSIRC    } from '../fli_tools/psirc'
include { JCCIRC   } from '../fli_tools/jccirc'
include { CIRCTOOLS } from '../fli_tools/circtools'

workflow FLI_DETECTION {
    take:
    ch_reads
    ch_reads_fixed_length
    ch_fasta
    ch_gtf
    ch_bwa_index
    ch_ciri2_txt
    ch_ciri2_sam
    ch_bsj_annotation
    ch_bsj_bed12_per_sample
    ch_bsj_reads
    ch_psirc_index
    ch_psirc_bsj
    ch_star_bam
    ch_star_junction
    ch_bsj_bed_per_sample
    fli_tools
    max_shift

    main:

    ch_bed12 = channel.empty()
    ch_fli_fasta = channel.empty()

    if (fli_tools.contains('cirifull')) {
        CIRIFULL(ch_reads_fixed_length, ch_bsj_annotation, ch_fasta, ch_gtf, ch_bwa_index, ch_ciri2_txt, ch_ciri2_sam, max_shift)
        ch_bed12 = ch_bed12.mix(CIRIFULL.out.bed12)
        ch_fli_fasta = ch_fli_fasta.mix(CIRIFULL.out.fasta)
    }

    if (fli_tools.contains('psirc')) {
        PSIRC(ch_reads, ch_psirc_bsj, ch_psirc_index)
        ch_bed12 = ch_bed12.mix(PSIRC.out.bed12)
        ch_fli_fasta = ch_fli_fasta.mix(PSIRC.out.fasta)
    }

    if (fli_tools.contains('jccirc')) {
        JCCIRC(ch_reads, ch_bsj_bed12_per_sample, ch_bsj_reads, ch_fasta, ch_gtf)
    }

    if (fli_tools.contains('circtools')) {
        CIRCTOOLS(ch_bsj_bed_per_sample, ch_star_bam, ch_star_junction, ch_fasta, ch_gtf)
        ch_bed12 = ch_bed12.mix(CIRCTOOLS.out.bed12)
        ch_fli_fasta = ch_fli_fasta.mix(CIRCTOOLS.out.fasta)
    }

    emit:
    bed12 = ch_bed12
    fasta = ch_fli_fasta
}
