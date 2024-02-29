include { HISAT2_ALIGN              } from '../../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_SORT             } from '../../modules/nf-core/samtools/sort/main'
include { STRINGTIE_STRINGTIE       } from '../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_PREPDE          } from '../../modules/local/stringtie/prepde/main'
include { PARENT_GENE               } from '../../modules/local/annotation/parent_gene/main'
include { PREPARE_CLR_TEST          } from '../../modules/local/circtest/prepare/main'
include { CIRCTEST                  } from '../../modules/local/circtest/test/main'
include { DESEQ2_DIFFERENTIAL_EXPRESSION   } from '../../modules/local/deseq2/differential_expression/main'
include { BAM_SORT_STATS_SAMTOOLS   } from '../nf-core/bam_sort_stats_samtools/main'

workflow DIFFERENTIAL_EXPRESSION {

    take:
    reads
    gtf
    fasta
    hisat2_index
    splice_sites
    phenotype
    dea_matrix
    clr_matrix
    species
    ensembl_map
    exon_boundary

    main:
    ch_versions = Channel.empty()
    qc_reports = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)
    ch_gtf   = Channel.fromPath(gtf)

    ch_fasta.map{ it -> [ [id: it.simpleName], [it] ]
    }.set{ fasta_tuple }

    ch_gtf.map{ it -> [ [id: it.simpleName], [it] ]
    }.set{ gtf_tuple }

    //
    // LINEAR RNA ALIGNEMT WORKFLOW:
    //

    HISAT2_ALIGN( reads, hisat2_index, splice_sites )
    BAM_SORT_STATS_SAMTOOLS( HISAT2_ALIGN.out.bam, fasta_tuple )
    STRINGTIE_STRINGTIE( BAM_SORT_STATS_SAMTOOLS.out.bam, gtf )
    STRINGTIE_PREPDE( STRINGTIE_STRINGTIE.out.transcript_gtf.map{ meta, gtf -> return [ gtf ] }.collect() )

    qc_reports = qc_reports.mix(HISAT2_ALIGN.out.summary.map{ meta, log -> log})
    qc_reports = qc_reports.mix(BAM_SORT_STATS_SAMTOOLS.out.stats.map{ meta, stats -> stats})
    qc_reports = qc_reports.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.map{ meta, flagstat -> flagstat})
    qc_reports = qc_reports.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.map{ meta, idxstats -> idxstats})


    //
    // Circular, Linear Differential Expression
    //

    DESEQ2_DIFFERENTIAL_EXPRESSION( STRINGTIE_PREPDE.out.gene_matrix, phenotype, dea_matrix, species, ensembl_map )

    //
    // CircRNA - Host Gene Ratio tests
    //

    ch_biotypes = Channel.fromPath("${projectDir}/bin/unwanted_biotypes.txt")

    PARENT_GENE( clr_matrix, gtf, ch_biotypes.collect(), exon_boundary )
    PREPARE_CLR_TEST( STRINGTIE_PREPDE.out.gene_matrix, clr_matrix, PARENT_GENE.out.circ_host_map, gtf )
    CIRCTEST( PREPARE_CLR_TEST.out.circular, PREPARE_CLR_TEST.out.linear, phenotype )

    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)
    ch_versions = ch_versions.mix(STRINGTIE_PREPDE.out.versions)
    ch_versions = ch_versions.mix(DESEQ2_DIFFERENTIAL_EXPRESSION.out.versions)
    ch_versions = ch_versions.mix(PARENT_GENE.out.versions)
    ch_versions = ch_versions.mix(PREPARE_CLR_TEST.out.versions)
    ch_versions = ch_versions.mix(CIRCTEST.out.versions)

    emit:
    versions = ch_versions
    reports = qc_reports
}
