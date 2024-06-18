include { BEDTOOLS_INTERSECT as INTERSECT           } from '../../../modules/nf-core/bedtools/intersect'
include { ANNOTATION as GTF_ANNOTATION              } from '../../../modules/local/annotation'
include { GAWK as INGEST_DATABASE_NAMES             } from '../../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_DATABASES             } from '../../../modules/nf-core/gnu/sort'
include { BEDTOOLS_INTERSECT as DATABASE_ANNOTATION } from '../../../modules/nf-core/bedtools/intersect'
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

    INTERSECT( regions.combine(gtf), [[], []] )
    GTF_ANNOTATION( INTERSECT.out.intersect, exon_boundary )

    ch_bed_merged = GTF_ANNOTATION.out.bed.filter{ meta, bed -> meta.tool == "merged" }
    ch_gtf_merged = GTF_ANNOTATION.out.gtf.filter{ meta, gtf -> meta.tool == "merged" }

    INGEST_DATABASE_NAMES( ch_annotation, [] )
    DATABASE_ANNOTATION( ch_bed_merged.combine(INGEST_DATABASE_NAMES.out.output)
        .map{ meta1, regions, meta2, database ->
            [meta1 + [id: "${meta1.id}:${meta2.id}",
                sample: meta1.id,
                db: meta2.id,
                min_overlap: meta2.min_overlap], regions, database] },
        [[], []])
    COMBINE_DATABASES( DATABASE_ANNOTATION.out.intersect
        .map{ meta, annotated -> [[id: meta.sample], annotated] }
        .groupTuple())

    COMBINE_BEDS(ch_bed_merged.map{ meta, bed -> bed}.collect().map{[[id: "annotation"], it]})
    REMOVE_SCORE_STRAND( COMBINE_BEDS.out.sorted, [])
    COMBINE_GTFS(ch_gtf_merged.map{ meta, gtf -> gtf}.collect().map{[[id: "annotation"], it]})

    ch_versions = ch_versions.mix(INTERSECT.out.versions)
    ch_versions = ch_versions.mix(GTF_ANNOTATION.out.versions)
    ch_versions = ch_versions.mix(INGEST_DATABASE_NAMES.out.versions)
    ch_versions = ch_versions.mix(COMBINE_DATABASES.out.versions)
    ch_versions = ch_versions.mix(DATABASE_ANNOTATION.out.versions)
    ch_versions = ch_versions.mix(COMBINE_BEDS.out.versions)
    ch_versions = ch_versions.mix(REMOVE_SCORE_STRAND.out.versions)
    ch_versions = ch_versions.mix(COMBINE_GTFS.out.versions)

    emit:
    merged_bed = ch_bed_merged
    bed = REMOVE_SCORE_STRAND.out.output
    gtf = COMBINE_GTFS.out.sorted

    versions = ch_versions
}
