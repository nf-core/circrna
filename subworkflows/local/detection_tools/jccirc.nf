include { TRINITY } from '../../../modules/nf-core/trinity'

workflow JCCIRC {
    take:
    reads

    main:
    ch_versions = Channel.empty()

    ch_trinity = reads.map{ _meta, r -> [[id: 'all_samples'], r] }.groupTuple()
    TRINITY(ch_trinity)
    // ch_versions = ch_versions.mix(TRINITY.out.versions)

    // JCCIRC(
    // )

    emit:
    versions = ch_versions
}
