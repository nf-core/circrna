include { TRINITY } from '../../../modules/nf-core/trinity'
include { JCCIRC_JCCIRC as MAIN } from '../../../modules/local/jccirc/jccirc'
include { JCCIRC_PREP as PREP } from '../../../modules/local/jccirc/prep'

workflow JCCIRC {
    take:
    reads
    ch_bsj_annotation
    ch_bsj_reads
    ch_fasta
    ch_gtf

    main:
    ch_versions = Channel.empty()

    ch_trinity = reads
        .map{ _meta, r -> [[id: 'all_samples'], r] }
        .groupTuple()
        .map{ meta, r -> [meta, r.flatten()] }
    TRINITY(ch_trinity)
    ch_versions = ch_versions.mix(TRINITY.out.versions)

    PREP(
        ch_bsj_annotation.join(ch_bsj_reads)
    )
    ch_versions = ch_versions.mix(PREP.out.versions)

    MAIN(
        reads
            .join(PREP.out.merged)
            .combine(TRINITY.out.transcript_fasta
                .map{ _meta, fasta -> fasta}
            ),
        ch_fasta,
        ch_gtf
    )

    emit:
    versions = ch_versions
}
