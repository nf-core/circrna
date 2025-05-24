include { CIRCTOOLS_ANNOTATION as ANNOTATION } from '../../../modules/local/circtools/annotation'
include { GAWK as FORMAT_READS               } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_SORT                      } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX                     } from '../../../modules/nf-core/samtools/index'
include { CIRCTOOLS_RECONSTRUCT              } from '../../../modules/local/circtools/reconstruct'

workflow CIRCTOOLS {
    take:
    ch_bsj_bed_per_sample
    ch_star_bam
    ch_star_junction
    ch_fasta
    ch_gtf

    main:
    ch_versions = Channel.empty()

    ANNOTATION(ch_gtf)
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)

    // FORMAT_READS(ch_bsj_reads, [], false)
    // ch_versions = ch_versions.mix(FORMAT_READS.out.versions)

    SAMTOOLS_SORT(ch_star_bam, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    CIRCTOOLS_RECONSTRUCT(
        ch_bsj_bed_per_sample.map { meta, bed ->  [[id: meta.id], bed] }
        .join(SAMTOOLS_SORT.out.bam.map { meta, bam ->  [[id: meta.id], bam] })
        .join(SAMTOOLS_INDEX.out.bai.map { meta, bai ->  [[id: meta.id], bai] })
        .join(ch_star_junction.map { meta, junction ->  [[id: meta.id], junction] }),
        ANNOTATION.out.bed,
    )
    ch_versions = ch_versions.mix(CIRCTOOLS_RECONSTRUCT.out.versions)

    emit:
    versions = ch_versions
}
