include { GNU_SORT as SORT      } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE        } from '../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_INTERSECT    } from '../../modules/nf-core/bedtools/intersect'
include { BEDTOOLS_JACCARD      } from '../../modules/nf-core/bedtools/jaccard'
include { BEDTOOLS_GENOMECOV    } from '../../modules/nf-core/bedtools/genomecov'
include { BENCHMARKING_MULTIQC  } from '../../modules/local/benchmarking/multiqc'
include { LOCATION_PLOT         } from '../../modules/local/benchmarking/location_plots'
include { OVERLAP_PLOT          } from '../../modules/local/benchmarking/overlap_plot'
include { SEQ_DEPTH_CORRELLATION} from '../../modules/local/benchmarking/seq_depth_plot'
include { WRITE                 } from '../../modules/local/benchmarking/write_file'
include { WRITE as WRITE2       } from '../../modules/local/benchmarking/write_file'


workflow BENCHMARKING {

    take:
    ch_real_bed
    ch_benchmarking_bed
    ch_real_bam
    ch_benchmarking_bam
    ch_trim_report

    main:

    ch_versions = Channel.empty()

    ch_all = ch_real_bed.mix(ch_benchmarking_bed)
        .map{ meta, bed -> [[id: meta.tool + "_" + (meta.benchmarking ? "benchmarking" : "real"),
                            tool: meta.tool,
                            benchmarking: meta.benchmarking], bed]}
        .groupTuple()

    ch_all.view {"all: $it"}

    SORT(ch_all)

    SORT.out.sorted.view { "sort: $it"}

    BEDTOOLS_MERGE(SORT.out.sorted).bed.branch{ meta, bed ->
            real: !meta.benchmarking
            benchmarking: meta.benchmarking
        }.set { ch_merged }

    ch_merged.real.view { "merged_rea: $it" }

    ch_merged.benchmarking.view { "merged_ben: $it" }

    ch_joined = ch_merged.real.map{ meta, bed -> [[id: meta.tool], bed]}
        .join(ch_merged.benchmarking.map{ meta, bed -> [[id: meta.tool], bed]})

    ch_joined.view { "joined1: $it"}
    ch_intersect = BEDTOOLS_INTERSECT(ch_joined,[[], []])

    OVERLAP_PLOT(ch_intersect)

    ch_joined.view { "joined2: $it"}

    LOCATION_PLOT(ch_joined)

    ch_joined.view { "joined3: $it"}

    ch_meta = ch_real_bam.map { it[0] }
    ch_path = ch_real_bam.map { it[1] }
    ch_scale = Channel.value(1)
    ch_genomecov_inputs = ch_meta.combine(ch_path).combine(ch_scale)
        .map { meta, path, scale ->
            tuple(meta, path, scale)
        }

    ch_genomecov = BEDTOOLS_GENOMECOV(ch_genomecov_inputs, [], "bg",false)

    ch_joined.view { "joined: $it"}

    ch_corr = SEQ_DEPTH_CORRELLATION(ch_joined, ch_genomecov.genomecov)
    ch_corr.view { "emits: $it"}


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

    ch_reports = BENCHMARKING_MULTIQC.out.report.mix(LOCATION_PLOT.out.report)
    ch_reports = ch_reports.mix(OVERLAP_PLOT.out.report)

    emit:
    reports = ch_reports
    versions       = ch_versions                     // channel: [ versions.yml ]
}