include { CIRIQUANT as RUN } from '../../../modules/local/ciriquant/ciriquant'
include { GAWK as UNIFY    } from '../../../modules/nf-core/gawk'

workflow CIRIQUANT {
    take:
    reads
    ch_gtf
    ch_fasta
    bwa_index
    hisat2_index
    
    main:
    ch_versions = Channel.empty()

    RUN( reads, ch_gtf, ch_fasta, bwa_index, hisat2_index )
    UNIFY( RUN.out.gtf.map{ meta, gtf ->
        [ meta + [tool: "ciriquant"], gtf ] }, [] )

    ch_versions = ch_versions.mix(RUN.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}