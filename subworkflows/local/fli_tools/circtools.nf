include { CIRCTOOLS_ANNOTATION as ANNOTATION } from '../../../modules/local/circtools/annotation'
include { GAWK as FORMAT_READS               } from '../../../modules/nf-core/gawk'
include { CIRCTOOLS_RECONSTRUCT              } from '../../../modules/local/circtools/reconstruct'

workflow CIRCTOOLS {
    take:
    ch_bsj_reads
    ch_star_bam
    ch_gtf

    main:
    ch_versions = Channel.empty()

    ANNOTATION(ch_gtf)
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)

    FORMAT_READS(ch_bsj_reads, [], false)
    ch_versions = ch_versions.mix(FORMAT_READS.out.versions)

    CIRCTOOLS_RECONSTRUCT(
        FORMAT_READS.out.output.map { meta, reads ->  [[id: meta.id], reads] }
        .join(ch_star_bam.map { meta, bam ->  [[id: meta.id], bam] }),
        ANNOTATION.out.bed,
    )
    ch_versions = ch_versions.mix(CIRCTOOLS_RECONSTRUCT.out.versions)

    emit:
    versions = ch_versions
}
