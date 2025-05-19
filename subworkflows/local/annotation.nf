include { GAWK as INGEST_DATABASE_NAMES            } from '../../modules/nf-core/gawk'
include { GNU_SORT as COMBINE_DATABASES            } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_INTERSECT as INTERSECT_DATABASE } from '../../modules/nf-core/bedtools/intersect'
include { CIRCEXPLORER2_ANNOTATE as ANNOTATE       } from '../../modules/nf-core/circexplorer2/annotate'
include { BEDTOOLS_GETFASTA as GET_FASTA           } from '../../modules/nf-core/bedtools/getfasta'
include { GAWK as RENAME                           } from '../../modules/nf-core/gawk'
include { GAWK as CUT_BED12                        } from '../../modules/nf-core/gawk'
include { ANNOTATION_BED2GTF as BED2GTF            } from '../../modules/local/annotation/bed2gtf'

workflow ANNOTATION {
    take:
    regions
    ch_annotation
    fasta
    circexplorer2_index

    main:
    ch_versions = Channel.empty()

    INGEST_DATABASE_NAMES(ch_annotation, [], false)
    ch_versions = ch_versions.mix(INGEST_DATABASE_NAMES.out.versions)

    INTERSECT_DATABASE(
        regions.combine(INGEST_DATABASE_NAMES.out.output).map { meta1, _regions, meta2, database ->
            [
                [
                    id: "${meta1.id}-${meta2.id}",
                    tool: meta1.tool,
                    original_meta: meta1,
                    min_overlap: meta2.min_overlap,
                ],
                _regions,
                database,
            ]
        },
        [[], []],
    )
    ch_versions = ch_versions.mix(INTERSECT_DATABASE.out.versions)

    ANNOTATE(regions, fasta, circexplorer2_index)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)

    RENAME(ANNOTATE.out.txt, [], false)
    ch_versions = ch_versions.mix(RENAME.out.versions)

    BED2GTF(RENAME.out.output.map { meta, bed12 -> [meta, bed12, []] }, params.exons_only)
    ch_versions = ch_versions.mix(BED2GTF.out.versions)

    CUT_BED12(RENAME.out.output, [], false)
    ch_versions = ch_versions.mix(CUT_BED12.out.versions)

    GET_FASTA(RENAME.out.output, fasta)
    ch_versions = ch_versions.mix(GET_FASTA.out.versions)

    emit:
    bed12    = RENAME.out.output
    gtf      = BED2GTF.out.gtf
    fasta    = GET_FASTA.out.fasta
    versions = ch_versions
}
