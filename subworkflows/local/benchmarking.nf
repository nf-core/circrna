include { GNU_SORT as SORT          } from '../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE            } from '../../modules/nf-core/bedtools/merge'
include { BEDTOOLS_INTERSECT        } from '../../modules/nf-core/bedtools/intersect'
include { BEDTOOLS_JACCARD          } from '../../modules/nf-core/bedtools/jaccard'
include { BEDTOOLS_GENOMECOV        } from '../../modules/nf-core/bedtools/genomecov'
include { BENCHMARKING_MULTIQC as JACCARD_MULTIQC     } from '../../modules/local/benchmarking/multiqc'
include { BENCHMARKING_MULTIQC as CORRELATION_MULTIQC } from '../../modules/local/benchmarking/multiqc'
include { PNG_JSON as LOCATION_JSON } from '../../modules/local/benchmarking/png_json'
include { PNG_JSON as OVERLAP_JSON  } from '../../modules/local/benchmarking/png_json'
include { LOCATION_PLOT             } from '../../modules/local/benchmarking/location_plots'
include { OVERLAP_PLOT              } from '../../modules/local/benchmarking/overlap_plot'
include { SEQ_DEPTH_CORRELLATION    } from '../../modules/local/benchmarking/seq_depth_plot'
include { AVERAGE_TSV               } from '../../modules/local/benchmarking/average_tsv'
include { POLYATAILOR as POLYATAILOR_REAL     } from '../../modules/local/benchmarking/polyatailor'
include { POLYATAILOR as POLYATAILOR_BENCHMARK} from '../../modules/local/benchmarking/polyatailor'
include { PLOT_POLYATAILS           } from '../../modules/local/benchmarking/plot_polyatails'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_TOTAL       } from '../../modules/nf-core/subread/featurecounts'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_BENCHMARKING} from '../../modules/nf-core/subread/featurecounts'
include { STAR2PASS as STAR2PASS_REAL            } from './detection_tools/star2pass'
include { STAR2PASS as STAR2PASS_BENCHMARKING    } from './detection_tools/star2pass'
include { STAR_GENOMEGENERATE       } from '../../modules/nf-core/star/genomegenerate'
include { DECOMPRESS_READS as DECOMPRESS_REAL        } from '../../modules/local/benchmarking/decompress_reads'
include { DECOMPRESS_READS as DECOMPRESS_BENCHMARKING} from '../../modules/local/benchmarking/decompress_reads'
include { RRNA_CORRELATION          } from '../../modules/local/benchmarking/rRNA_correlation'



workflow BENCHMARKING {

    take:
    ch_reads_real
    ch_reads_benchmarking
    ch_real_bed
    ch_benchmarking_bed
    ch_real_bam
    ch_benchmarking_bam
    ch_trim_report
    ch_fasta
    bsj_reads

    main:

    //quality score estimation
    /*ch_reads_benchmarking.view {"rb: $it"}
    ch_reads_real.view {"rr: $it"}
    ch_polya_real = POLYATAILOR_REAL(ch_reads_real, ch_real_bam.map { it[1] })
    ch_polya_benchmarking = POLYATAILOR_BENCHMARK(ch_reads_benchmarking, ch_benchmarking_bam.map { it[1] })

    ch_polya_real.tails.collect().view {"pr: $it"}
    ch_polya_benchmarking.tails.collect().view {"pb: $it"}

    ch_polya_plots = PLOT_POLYATAILS(ch_polya_real.tails.collect(), ch_polya_benchmarking.tails.collect())
    ch_polya_plots.average_tails.view {"$it"}*/


    // Use only the paths in each tuple of ch_reads_real and ch_reads_benchmarking, assuming the metadata is already present
    // Define parameters
    ch_benchmark_gtf = "/nfs/data3/CIRCEST/runs/test_benchmarking/gencode.v47.primary_assembly.annotation.gtf"
    ch_benchmarking_fasta = "/nfs/data3/CIRCEST/runs/test_benchmarking/gencode.v47.rRNARNA_transcripts.fa"
    ch_filter_gtf = "/nfs/data3/CIRCEST/runs/test_benchmarking/ensembl_rRNA.gtf"

    STAR_GENOMEGENERATE(ch_fasta,  tuple([id: "benchmarking_gtf"], file(ch_benchmark_gtf)))
    star_index          = params.star     ? Channel.value([[id: "star"], file(params.star, checkIfExists: true)])       : STAR_GENOMEGENERATE.out.index.collect()
    star_ignore_sjdbgtf = true
    seq_center = params.seq_center ?: ''
    seq_platform = ''


    ch_reads_real_restructured = ch_reads_real.map { meta, fastq_gz_list ->
    tuple(meta, fastq_gz_list[0], fastq_gz_list[1])
    }
    ch_reads_benchmarking_restructured = ch_reads_benchmarking.map { meta, fastq_gz_list ->
    tuple(meta, fastq_gz_list[0], fastq_gz_list[1])
    }
    ch_uncompressed_reads_real = DECOMPRESS_REAL(ch_reads_real_restructured)
    ch_uncompressed_reads_benchmarking = DECOMPRESS_BENCHMARKING(ch_reads_benchmarking_restructured)
    
    ch_rRNA_real_bam = STAR2PASS_REAL(ch_uncompressed_reads_real, star_index, tuple([id: "Benchmarking_gtf"], file(ch_benchmark_gtf)), bsj_reads, star_ignore_sjdbgtf, seq_center, seq_platform).bam

    ch_rRNA_benchmarking_bam = STAR2PASS_BENCHMARKING(ch_uncompressed_reads_benchmarking, star_index, tuple([id: "Benchmarking_gtf"], file(ch_benchmark_gtf)), bsj_reads, star_ignore_sjdbgtf, seq_center, seq_platform).bam



    ch_rRNA_real_input = ch_real_bam.map { meta, path ->
    tuple(meta, path, file(ch_filter_gtf))
    }
    ch_rRNA_benchmarking_input = ch_benchmarking_bam.map { meta, path ->
    tuple(meta, path, file(ch_filter_gtf))
    }
    
    ch_rRNA_real = FEATURECOUNTS_TOTAL(ch_rRNA_real_input).summary
    ch_rRNA_benchmarking = FEATURECOUNTS_BENCHMARKING(ch_rRNA_benchmarking_input).summary


    ch_collected_real_bed = ch_real_bed.collect()
    ch_collected_benchmarking_bed = ch_benchmarking_bed.collect()
    ch_collected_rRNA_real = ch_rRNA_real.collect()
    ch_collected_rRNA_benchmarking = ch_rRNA_benchmarking.collect()

    RRNA_CORRELATION(
        ch_collected_real_bed,
        ch_collected_benchmarking_bed,
        ch_collected_rRNA_real,
        ch_collected_rRNA_benchmarking
    )


    //ch_rRNA_benchmarking.view { "bench: $it" }


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
    OVERLAP_PLOT(ch_intersect.intersect)
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
