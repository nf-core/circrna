include { CIRCTOOLS_ANNOTATION as ANNOTATION   } from '../../../modules/local/circtools/annotation'
include { SAMTOOLS_SORT                        } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX                       } from '../../../modules/nf-core/samtools/index'
include { CIRCTOOLS_RECONSTRUCT as RECONSTRUCT } from '../../../modules/local/circtools/reconstruct'
include { BEDTOOLS_GETFASTA                    } from '../../../modules/nf-core/bedtools/getfasta'

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

    SAMTOOLS_SORT(ch_star_bam, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    RECONSTRUCT(
        ch_bsj_bed_per_sample.map { meta, bed -> [[id: meta.id], bed] }
        .join(SAMTOOLS_SORT.out.bam.map { meta, bam -> [[id: meta.id], bam] })
        .join(SAMTOOLS_INDEX.out.bai.map { meta, bai -> [[id: meta.id], bai] })
        .join(ch_star_junction.map { meta, junction -> [[id: meta.id], junction] }),
        ANNOTATION.out.bed,
    )
    ch_versions = ch_versions.mix(RECONSTRUCT.out.versions)

    BEDTOOLS_GETFASTA(RECONSTRUCT.out.exon_counts_bed, ch_fasta.map { _meta, fasta -> fasta })
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

    emit:
    fasta = BEDTOOLS_GETFASTA.out.fasta.map{ meta, fasta -> [meta + [fli_tool: 'circtools'], fasta] }
    bed12 = RECONSTRUCT.out.exon_counts_bed.map{ meta, bed -> [meta + [fli_tool: 'circtools'], bed] }

    versions = ch_versions
}
