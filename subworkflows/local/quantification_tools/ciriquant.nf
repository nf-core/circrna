include { CIRIQUANT as MAIN                } from '../../../modules/local/ciriquant'
include { PYGTFTK_TABULATE as EXTRACT_CIRC } from '../../../modules/local/pygtftk/tabulate'
include { GAWK as EXTRACT_GENES            } from '../../../modules/nf-core/gawk'
include { GAWK as NAME_EXPRESSION          } from '../../../modules/nf-core/gawk'
include { CSVTK_JOIN as JOIN_GENE          } from '../../../modules/nf-core/csvtk/join'
include { CSVTK_JOIN as JOIN_CIRC          } from '../../../modules/nf-core/csvtk/join'


workflow CIRIQUANT {
    take:
    reads
    ch_bed
    ch_gtf
    ch_fasta
    bwa_index
    hisat2_index

    main:
    ch_versions = Channel.empty()

    MAIN( reads, ch_bed, ch_gtf, ch_fasta, bwa_index, hisat2_index )
    ch_versions = ch_versions.mix(MAIN.out.versions)

    EXTRACT_CIRC( MAIN.out.gtf )
    ch_versions = ch_versions.mix(EXTRACT_CIRC.out.versions)

    EXTRACT_GENES( MAIN.out.gene, [] )
    ch_versions = ch_versions.mix(EXTRACT_GENES.out.versions)

    ch_tables = EXTRACT_CIRC.out.table
        .map{ meta, table -> [meta + [type: 'circ'], table] }
        .mix(EXTRACT_GENES.out.output
            .map{ meta, table -> [meta + [type: 'gene'], table] }
        )
    
    NAME_EXPRESSION(ch_tables, [])
    ch_versions = ch_versions.mix(NAME_EXPRESSION.out.versions)

    ch_tables = NAME_EXPRESSION.out.output.branch{ meta, table ->
        circ: meta.type == 'circ'
        gene: meta.type == 'gene'    
    }

    JOIN_GENE(
        ch_tables.gene.map{meta, table -> table}.collect().map{[[id: "gene"], it]}
    )
    ch_versions = ch_versions.mix(JOIN_GENE.out.versions)

    JOIN_CIRC(
        ch_tables.circ.map{meta, table -> table}.collect().map{[[id: "circ"], it]}
    )
    ch_versions = ch_versions.mix(JOIN_CIRC.out.versions)
    
    emit:
    gene_tpm = JOIN_GENE.out.csv
    circ_cpm = JOIN_CIRC.out.csv

    versions = ch_versions
}
