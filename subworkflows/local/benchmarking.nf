include { GNU_SORT as SORT } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE   } from '../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_JACCARD } from '../../modules/nf-core/bedtools/jaccard'

workflow BENCHMARKING {

    take:
    ch_real_bed
    ch_benchmarking_bed

    main:

    ch_versions = Channel.empty()

    ch_all = ch_real_bed.mix(ch_benchmarking_bed)
        .map{ meta, bed -> [[id: meta.tool + "_" + (meta.benchmarking ? "benchmarking" : "real"),
                            tool: meta.tool,
                            benchmarking: meta.benchmarking], bed]}
        .groupTuple()
    
    SORT(ch_all)
    BEDTOOLS_MERGE(SORT.out.sorted).bed.branch{ meta, bed -> 
            real: !meta.benchmarking
            benchmarking: meta.benchmarking
        }.set { ch_merged }

    ch_joined = ch_merged.real.map{ meta, bed -> [[id: meta.tool], bed]}
        .join(ch_merged.benchmarking.map{ meta, bed -> [[id: meta.tool], bed]})

    ch_jaccard = BEDTOOLS_JACCARD(ch_joined, [[], []]).tsv

    ch_versions = ch_versions.mix(SORT.out.versions)
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)
    ch_versions = ch_versions.mix(BEDTOOLS_JACCARD.out.versions)

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
