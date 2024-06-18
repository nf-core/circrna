include { BEDTOOLS_INTERSECT as INTERSECT_GTF       } from '../../../modules/nf-core/bedtools/intersect'
include { ANNOTATION as ANNOTATE                    } from '../../../modules/local/annotation'
include { GAWK as INGEST_DATABASE_NAMES             } from '../../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_DATABASES             } from '../../../modules/nf-core/gnu/sort'
include { BEDTOOLS_INTERSECT as INTERSECT_DATABASE  } from '../../../modules/nf-core/bedtools/intersect'
include { GAWK as CLEAN_DATABASE_ANNOTATION         } from '../../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_BEDS                  } from '../../../modules/nf-core/gnu/sort'
include { GAWK as REMOVE_SCORE_STRAND               } from '../../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_GTFS                  } from '../../../modules/nf-core/gnu/sort'

workflow ANNOTATION {
    take:
    regions
    gtf
    exon_boundary
    ch_annotation

    main:
    ch_versions = Channel.empty()

    INTERSECT_GTF( regions.combine(gtf), [[], []] )
    INGEST_DATABASE_NAMES( ch_annotation, [] )
    INTERSECT_DATABASE( regions.combine(INGEST_DATABASE_NAMES.out.output)
        .map{ meta1, regions, meta2, database ->
            [[id: "${meta1.id}-${meta2.id}",
                original_meta: meta1,
                min_overlap: meta2.min_overlap], regions, database] },
        [[], []])

    ANNOTATE( INTERSECT_GTF.out.intersect
        .join(INTERSECT_DATABASE.out.intersect
            .map{ meta, bed -> [meta.original_meta, bed] }
            .groupTuple(), remainder: true)
        .map{ meta, gtf_intersection, db_intersections -> [meta, gtf_intersection, db_intersections ?: []]},
        exon_boundary )

    ch_bed_merged = ANNOTATE.out.bed.filter{ meta, bed -> meta.tool == "merged" }
    ch_gtf_merged = ANNOTATE.out.gtf.filter{ meta, gtf -> meta.tool == "merged" }

    COMBINE_BEDS(ch_bed_merged.map{ meta, bed -> bed}.collect().map{[[id: "annotation"], it]})
    REMOVE_SCORE_STRAND( COMBINE_BEDS.out.sorted, [])
    COMBINE_GTFS(ch_gtf_merged.map{ meta, gtf -> gtf}.collect().map{[[id: "annotation"], it]})

    ch_versions = ch_versions.mix(INTERSECT_GTF.out.versions)
    ch_versions = ch_versions.mix(INGEST_DATABASE_NAMES.out.versions)
    ch_versions = ch_versions.mix(INTERSECT_DATABASE.out.versions)
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
