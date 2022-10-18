

include { STAR_GENOMEGENERATE } from '../../modules/nf-core/modules/star/genomegenerate/main'

workflow PREPARE_GENOME {

    take:
    fasta
    gtf

    main:
    ch_versions = Channel.empty()

    STAR_GENOMEGENERATE(
        fasta,
        gtf
    )

    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    star     = STAR_GENOMEGENERATE.out.index
    versions = ch_versions
}
