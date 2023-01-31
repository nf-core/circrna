include { HISAT2_ALIGN              } from '../../modules/nf-core/hisat2/align/main'
include { SAMTOOLS_SORT             } from '../../modules/nf-core/samtools/sort/main'
include { STRINGTIE_STRINGTIE       } from '../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_PREPDE          } from '../../modules/local/stringtie/prepde/main'
include { DESEQ2_DIFFERENTIAL_EXPRESSION   } from '../../modules/local/deseq2/differential_expression/main'
include { PARENT_GENE               } from '../../modules/local/annotation/parent_gene/main'
include { PREPARE_CLR_TEST          } from '../../modules/local/circtest/prepare/main'
include { CIRCTEST                  } from '../../modules/local/circtest/test/main'

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

    //
    // LINEAR RNA ALIGNEMT WORKFLOW:
    //

    HISAT2_ALIGN( reads, hisat2_index, splice_sites )
    SAMTOOLS_SORT( HISAT2_ALIGN.out.bam )
    STRINGTIE_STRINGTIE( SAMTOOLS_SORT.out.bam, gtf )
    STRINGTIE_PREPDE( STRINGTIE_STRINGTIE.out.transcript_gtf.map{ meta, gtf -> return [ gtf ] }.collect() )

    //
    // Circular, Linear Differential Expression
    //

    DESEQ2_DIFFERENTIAL_EXPRESSION( STRINGTIE_PREPDE.out.gene_matrix, phenotype, dea_matrix, species, ensembl_map )

    //
    // CircRNA - Host Gene Ratio tests
    //

    PARENT_GENE( clr_matrix, gtf, exon_boundary )
    PREPARE_CLR_TEST( STRINGTIE_PREPDE.out.gene_matrix, clr_matrix, PARENT_GENE.out.circ_host_map, gtf )
    CIRCTEST( PREPARE_CLR_TEST.out.circular, PREPARE_CLR_TEST.out.linear, phenotype )

    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)
    ch_versions = ch_versions.mix(STRINGTIE_PREPDE.out.versions)
    ch_versions = ch_versions.mix(DESEQ2_DIFFERENTIAL_EXPRESSION.out.versions)
    ch_versions = ch_versions.mix(PARENT_GENE.out.versions)
    ch_versions = ch_versions.mix(PREPARE_CLR_TEST.out.versions)
    ch_versions = ch_versions.mix(CIRCTEST.out.versions)

    emit:
    versions = ch_versions
}
