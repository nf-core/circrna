include { CIRIQUANT as RUN           } from '../../../modules/local/ciriquant/ciriquant'
include { CIRIQUANT_FILTER as FILTER } from '../../../modules/local/ciriquant/filter'

workflow CIRIQUANT {
    take:
    reads
    ch_gtf
    ch_fasta
    bwa_index
    hisat2_index
    bsj_reads
    
    main:
    ch_versions = Channel.empty()

    RUN( reads, ch_gtf, ch_fasta, bwa_index, hisat2_index )
    FILTER( RUN.out.gtf.map{ meta, gtf ->
        [ meta + [tool: "ciriquant"], gtf ] }, bsj_reads )

    ch_versions = ch_versions.mix(RUN.out.versions)
    ch_versions = ch_versions.mix(FILTER.out.versions)

    emit:
    matrix = FILTER.out.matrix
    results = FILTER.out.results

    versions = ch_versions
}