include { CIRIFULL } from './fli_tools/cirifull'
include { PSIRC    } from './fli_tools/psirc'
include { JCCIRC   } from './fli_tools/jccirc'

workflow FLI_DETECTION {
    take:
    ch_reads
    ch_reads_fixed_length
    ch_fasta
    ch_gtf
    ch_bwa_index
    ch_ciri_txt
    ch_ciri_sam
    ch_bsj_annotation
    ch_bsj_reads
    ch_psirc_index
    ch_psirc_bsj

    main:
    ch_versions = Channel.empty()

    def fli_tools = params.fli_tools.split(',').collect { it.trim() }

    if (fli_tools.contains('cirifull')) {
        CIRIFULL(ch_reads_fixed_length, ch_bsj_annotation, ch_fasta, ch_gtf, ch_bwa_index, ch_ciri_txt, ch_ciri_sam)
        ch_versions = ch_versions.mix(CIRIFULL.out.versions)
    }

    if (fli_tools.contains('psirc')) {
        PSIRC(ch_reads, ch_psirc_bsj, ch_psirc_index)
        ch_versions = ch_versions.mix(PSIRC.out.versions)
    }

    if (fli_tools.contains('jccirc')) {
        JCCIRC(ch_reads, ch_bsj_annotation, ch_bsj_reads, ch_fasta, ch_gtf)
        ch_versions = ch_versions.mix(JCCIRC.out.versions)
    }

    emit:
    versions = ch_versions
}
