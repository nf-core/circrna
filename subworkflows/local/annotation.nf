include { BEDTOOLS_INTERSECT as INTERSECT_GTF       } from '../../modules/nf-core/bedtools/intersect'
include { GAWK as INGEST_DATABASE_NAMES             } from '../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_DATABASES             } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_INTERSECT as INTERSECT_DATABASE  } from '../../modules/nf-core/bedtools/intersect'
include { ANNOTATION as ANNOTATE                    } from '../../modules/local/annotation'

workflow ANNOTATION {
    take:
    regions
    ch_gtf
    ch_annotation

    main:
    ch_versions = Channel.empty()

    INTERSECT_GTF( regions.combine(ch_gtf.map{meta, gtf -> gtf}), [[], []] )
    ch_versions = ch_versions.mix(INTERSECT_GTF.out.versions)

    INGEST_DATABASE_NAMES( ch_annotation, [] )
    ch_versions = ch_versions.mix(INGEST_DATABASE_NAMES.out.versions)

    INTERSECT_DATABASE( regions.combine(INGEST_DATABASE_NAMES.out.output)
        .map{ meta1, regions, meta2, database ->
            [[id: "${meta1.id}-${meta2.id}",
                tool: meta1.tool,
                original_meta: meta1,
                min_overlap: meta2.min_overlap], regions, database] },
        [[], []])
    ch_versions = ch_versions.mix(INTERSECT_DATABASE.out.versions)

    ANNOTATE( INTERSECT_GTF.out.intersect
        .join(INTERSECT_DATABASE.out.intersect
            .map{ meta, bed -> [meta.original_meta, bed] }
            .groupTuple(), remainder: true)
        .map{ meta, gtf_intersection, db_intersections -> [meta, gtf_intersection, db_intersections ?: []]})
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)

    emit:
    bed = ANNOTATE.out.bed
    gtf = ANNOTATE.out.gtf

    versions = ch_versions
}
