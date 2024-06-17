include { BEDTOOLS_INTERSECT as INTERSECT } from '../../../modules/nf-core/bedtools/intersect'
include { ANNOTATION as ANNOTATE          } from '../../../modules/local/annotation'
include { GNU_SORT as COMBINE_BEDS        } from '../../../modules/nf-core/gnu/sort'
include { GAWK as REMOVE_SCORE_STRAND     } from '../../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_GTFS        } from '../../../modules/nf-core/gnu/sort'

workflow ANNOTATION {
    take:
    regions
    gtf
    exon_boundary
    ch_annotations

    main:
    ch_versions = Channel.empty()

    INTERSECT( regions.combine(gtf), [[], []])
    ANNOTATE( INTERSECT.out.intersect, exon_boundary )

    ch_annotations.view()

    ch_bed_merged = ANNOTATE.out.bed.filter{ meta, bed -> meta.tool == "merged" }
    ch_gtf_merged = ANNOTATE.out.gtf.filter{ meta, gtf -> meta.tool == "merged" }

    COMBINE_BEDS(ch_bed_merged.map{ meta, bed -> bed}.collect().map{[[id: "annotation"], it]})
    REMOVE_SCORE_STRAND( COMBINE_BEDS.out.sorted, [])
    COMBINE_GTFS(ch_gtf_merged.map{ meta, gtf -> gtf}.collect().map{[[id: "annotation"], it]})

    ch_versions = ch_versions.mix(INTERSECT.out.versions)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(COMBINE_BEDS.out.versions)
    ch_versions = ch_versions.mix(REMOVE_SCORE_STRAND.out.versions)
    ch_versions = ch_versions.mix(COMBINE_GTFS.out.versions)

    emit:
    merged_bed = ch_bed_merged
    bed = REMOVE_SCORE_STRAND.out.output
    gtf = COMBINE_GTFS.out.sorted

    versions = ch_versions
}
