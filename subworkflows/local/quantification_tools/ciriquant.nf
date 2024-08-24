include { CIRIQUANT as MAIN                } from '../../../modules/local/ciriquant'
include { PYGTFTK_TABULATE as EXTRACT_CIRC } from '../../../modules/local/pygtftk/tabulate'
include { GAWK as EXTRACT_GENES            } from '../../../modules/nf-core/gawk'
include { GAWK as NAME_EXPRESSION           } from '../../../modules/nf-core/gawk'

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
    
    emit:


    versions = ch_versions
}
