include { GNU_SORT as SORT     } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE       } from '../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_JACCARD     } from '../../modules/nf-core/bedtools/jaccard'
include { BENCHMARKING_MULTIQC } from '../../modules/local/benchmarking/multiqc'
include { PLOT_LOCI            } from '../../modules/local/plot_loci'


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

    PLOT_LOCI(ch_merged.real)

    ch_joined = ch_merged.real.map{ meta, bed -> [[id: meta.tool], bed]}
        .join(ch_merged.benchmarking.map{ meta, bed -> [[id: meta.tool], bed]})

    ch_jaccard = BEDTOOLS_JACCARD(ch_joined, [[], []]).tsv

    ch_stats = ch_jaccard.splitCsv(header: true, sep: "\t")
        .map{ meta, values -> [meta.id, values.intersection, values.union, values.jaccard, values.n_intersections]}
        .collectFile( newLine: true,
                        storeDir: params.outdir,
                        seed: "tool\tintersection\tunion\tjaccard\tn_intersections") {
                            row -> ["jaccard.tsv", row.join("\t")]
        }

    BENCHMARKING_MULTIQC(ch_stats)

    ch_versions = ch_versions.mix(SORT.out.versions)
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)
    ch_versions = ch_versions.mix(BEDTOOLS_JACCARD.out.versions)
    ch_versions = ch_versions.mix(BENCHMARKING_MULTIQC.out.versions)

    emit:
    reports = BENCHMARKING_MULTIQC.out.report
    versions       = ch_versions                     // channel: [ versions.yml ]
}
