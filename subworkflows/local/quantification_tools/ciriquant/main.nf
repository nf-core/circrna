include { CIRIQUANT as MAIN                } from '../../../../modules/local/ciriquant/ciriquant'
include { PYGTFTK_TABULATE as EXTRACT_CIRC } from '../../../../modules/local/pygtftk/tabulate'
include { GAWK as EXTRACT_GENES            } from '../../../../modules/nf-core/gawk'
include { JOIN_SAMPLES as JOIN_GENE        } from '../../../../modules/local/matrix/join_samples'
include { JOIN_SAMPLES as JOIN_CIRC        } from '../../../../modules/local/matrix/join_samples'

workflow CIRIQUANT {
    take:
    reads
    ch_bed
    ch_gtf
    ch_fasta
    bwa_index
    hisat2_index

    main:

    MAIN( reads, ch_bed, ch_gtf, ch_fasta, bwa_index, hisat2_index )

    EXTRACT_CIRC( MAIN.out.gtf )

    EXTRACT_GENES( MAIN.out.gene_list.map { meta, gene_list -> [meta, gene_list, 'gene.tsv'] }, [], false )

    JOIN_GENE(
        EXTRACT_GENES.out.output.map{meta, table -> [[id: 'gene'], meta.id, table]}.groupTuple(),
        "gene_id,gene_name"
    )

    JOIN_CIRC(
        EXTRACT_CIRC.out.table.map{meta, table -> [[id: 'circ'], meta.id, table]}.groupTuple(),
        "circ_id,gene_id"
    )

    emit:
    gene_tpm  = JOIN_GENE.out.joined
    circ_cpm  = JOIN_CIRC.out.joined
    raw       = MAIN.out.gtf
    stringtie = MAIN.out.gene_gtf

}
