include { TRINITY } from '../../../modules/nf-core/trinity'

workflow JCCIRC {
    take:
    reads
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

    JCCIRC(
        reads
            .join(ch_bsj_reads)
            .combine(TRINITY.out.transcript_fasta
                .map{ _meta, fasta -> fasta}
            ),
        ch_fasta,
        ch_gtf
    )

    emit:
    versions = ch_versions
}
