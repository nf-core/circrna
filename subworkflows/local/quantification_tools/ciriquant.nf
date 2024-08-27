include { CIRIQUANT as MAIN                } from '../../../modules/local/ciriquant/ciriquant'
include { PYGTFTK_TABULATE as EXTRACT_CIRC } from '../../../modules/local/pygtftk/tabulate'
include { GAWK as EXTRACT_GENES            } from '../../../modules/nf-core/gawk'
include { JOIN_SAMPLES as JOIN_GENE        } from '../../../modules/local/matrix/join_samples'
include { JOIN_SAMPLES as JOIN_CIRC        } from '../../../modules/local/matrix/join_samples'

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

    EXTRACT_GENES( MAIN.out.gene_list, [] )
    ch_versions = ch_versions.mix(EXTRACT_GENES.out.versions)

    JOIN_GENE(
        EXTRACT_GENES.out.output.map{meta, table -> [[id: 'gene'], meta.id, table]}.groupTuple()
    )
    ch_versions = ch_versions.mix(JOIN_GENE.out.versions)

    JOIN_CIRC(
        EXTRACT_CIRC.out.table.map{meta, table -> [[id: 'circ'], meta.id, table]}.groupTuple()
    )
    ch_versions = ch_versions.mix(JOIN_CIRC.out.versions)
    
    emit:
    gene_tpm  = JOIN_GENE.out.joined
    circ_cpm  = JOIN_CIRC.out.joined
    raw       = MAIN.out.gtf
    stringtie = MAIN.out.gene_gtf

    versions  = ch_versions
}
