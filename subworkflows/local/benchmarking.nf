include { GNU_SORT as SORT          } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE            } from '../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_INTERSECT        } from '../../modules/nf-core/bedtools/intersect'
include { BEDTOOLS_JACCARD          } from '../../modules/nf-core/bedtools/jaccard'
include { BEDTOOLS_GENOMECOV        } from '../../modules/nf-core/bedtools/genomecov'
include { BENCHMARKING_MULTIQC as JACCARD_MULTIQC } from '../../modules/local/benchmarking/multiqc'
include { BENCHMARKING_MULTIQC as CORRELATION_MULTIQC } from '../../modules/local/benchmarking/multiqc'
include { PNG_JSON as LOCATION_JSON } from '../../modules/local/benchmarking/png_json'
include { PNG_JSON as OVERLAP_JSON  } from '../../modules/local/benchmarking/png_json'
include { LOCATION_PLOT             } from '../../modules/local/benchmarking/location_plots'
include { OVERLAP_PLOT              } from '../../modules/local/benchmarking/overlap_plot'
include { SEQ_DEPTH_CORRELLATION    } from '../../modules/local/benchmarking/seq_depth_plot'
include { AVERAGE_TSV              } from '../../modules/local/benchmarking/average_tsv'



workflow BENCHMARKING {

    take:
    ch_real_bed
    ch_benchmarking_bed
    ch_real_bam
    ch_benchmarking_bam
    ch_trim_report

    main:

    //data preparation
    ch_versions = Channel.empty()

    ch_all = ch_real_bed.mix(ch_benchmarking_bed)
        .map{ meta, bed -> [[id: meta.tool + "_" + (meta.benchmarking ? "benchmarking" : "real"),
                            tool: meta.tool,
                            benchmarking: meta.benchmarking], bed]}
        .groupTuple()

    SORT(ch_all)
    ch_versions = ch_versions.mix(SORT.out.versions)

    BEDTOOLS_MERGE(SORT.out.sorted).bed.branch{ meta, bed ->
            real: !meta.benchmarking
            benchmarking: meta.benchmarking
        }.set { ch_merged }
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

    ch_joined = ch_merged.real.map{ meta, bed -> [[id: meta.tool], bed]}
        .join(ch_merged.benchmarking.map{ meta, bed -> [[id: meta.tool], bed]})

    //Overlap plot
    ch_intersect = BEDTOOLS_INTERSECT(ch_joined,[[], []])
    OVERLAP_PLOT(ch_intersect.out.intersect)
    OVERLAP_JSON(OVERLAP_PLOT.out.plots, "Overlap plots", "Plot the overlap circRNAs found in total and polyA data for the tools")
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)
    ch_versions = ch_versions.mix(OVERLAP_PLOT.out.versions)
    ch_versions = ch_versions.mix(OVERLAP_JSON.out.versions)

    //Location plot workflow
    LOCATION_PLOT(ch_joined)
    LOCATION_JSON(LOCATION_PLOT.out.plots, "Location plots", "Plots the location of the circRNAs found" )
    ch_versions = ch_versions.mix(LOCATION_PLOT.out.versions)
    ch_versions = ch_versions.mix(LOCATION_JSON.out.versions)

    //Pearson correllation workflow
    ch_meta = ch_real_bam.map { it[0] }
    ch_path = ch_real_bam.map { it[1] }
    ch_scale = Channel.value(1)
    ch_genomecov_inputs = ch_meta.combine(ch_path).combine(ch_scale)
        .map { meta, path, scale ->
            tuple(meta, path, scale)
        }

    ch_genomecov = BEDTOOLS_GENOMECOV(ch_genomecov_inputs, [], "bg",false)
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    ch_seqdepths = ch_genomecov.genomecov
        .map { genomecov_result -> genomecov_result[1].toString() }
        .collectFile(name: 'genomecov_paths.txt',
                        newLine: true)

    ch_corr = SEQ_DEPTH_CORRELLATION(ch_real_bed, ch_seqdepths.collect()).report
    ch_versions = ch_versions.mix(SEQ_DEPTH_CORRELLATION.out.versions)

    ch_pearson = ch_corr.splitCsv(header: true, sep: "\t")
    .map{ values -> [values.tool, values.pearson_corr]}
    .collectFile( newLine: true,
                    storeDir: params.outdir,
                    seed: "tool\tpearson_corr") {
                        row -> ["pearson.tsv", row.join("\t")]
    }

    AVERAGE_TSV(ch_pearson)
    CORRELATION_MULTIQC(AVERAGE_TSV.out.tsv)
    ch_versions = ch_versions.mix(AVERAGE_TSV.out.versions)
    ch_versions = ch_versions.mix(CORRELATION_MULTIQC.out.versions)

    //Jaccard Workflow
    ch_jaccard = BEDTOOLS_JACCARD(ch_joined, [[], []]).tsv
    ch_versions = ch_versions.mix(BEDTOOLS_JACCARD.out.versions)

    ch_stats = ch_jaccard.splitCsv(header: true, sep: "\t")
        .map{ meta, values -> [meta.id, values.intersection, values.union, values.jaccard, values.n_intersections]}
        .collectFile( newLine: true,
                        storeDir: params.outdir,
                        seed: "tool\tintersection\tunion\tjaccard\tn_intersections") {
                            row -> ["jaccard.tsv", row.join("\t")]
        }

    JACCARD_MULTIQC(ch_stats)
    ch_versions = ch_versions.mix(JACCARD_MULTIQC.out.versions)

    //combine results














    ch_reports = JACCARD_MULTIQC.out.report.mix(LOCATION_JSON.out.report)
    ch_reports = ch_reports.mix(OVERLAP_JSON.out.report)
    ch_reports = ch_reports.mix(CORRELATION_MULTIQC.out.report)

    emit:
    reports = ch_reports
    versions       = ch_versions                     // channel: [ versions.yml ]
}
