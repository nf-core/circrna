include { STAR2PASS as MATE1_STAR2PASS } from './star2pass'
include { STAR2PASS as MATE2_STAR2PASS } from './star2pass'
include { CIRCTOOLS_DETECT as DETECT   } from '../../../modules/local/circtools/detect'
include { UNIFY                        } from '../../../modules/local/circtools/unify'

workflow CIRCTOOLS {
    take:
    reads
    ch_fasta
    ch_gtf
    star_index
    star_junction
    ignore_sjdbgtf
    seq_platform
    seq_center
    bsj_reads

    main:
    ch_versions = Channel.empty()

    // Process mate 1
    ch_mate1 = reads
        .filter { meta, _reads -> !meta.single_end }
        .map { meta, _reads ->
            return [[id: meta.id, single_end: true], _reads[0]]
        }

    MATE1_STAR2PASS(
        ch_mate1,
        star_index,
        ch_gtf,
        bsj_reads,
        ignore_sjdbgtf,
        seq_center,
        seq_platform,
    )
    ch_versions = ch_versions.mix(MATE1_STAR2PASS.out.versions)

    // Process mate 2
    ch_mate2 = reads
        .filter { meta, _reads -> !meta.single_end }
        .map { meta, _reads ->
            return [[id: meta.id, single_end: true], _reads[1]]
        }

    MATE2_STAR2PASS(
        ch_mate2,
        star_index,
        ch_gtf,
        bsj_reads,
        ignore_sjdbgtf,
        seq_center,
        seq_platform,
    )
    ch_versions = ch_versions.mix(MATE2_STAR2PASS.out.versions)

    ch_combined_junctions = star_junction
        .map { meta, junction ->
            return [meta.id, meta, junction]
        }
        .join(
            MATE1_STAR2PASS.out.junction.map { meta, junction ->
                return [meta.id, junction]
            },
            remainder: true
        )
        .join(
            MATE2_STAR2PASS.out.junction.map { meta, junction ->
                return [meta.id, junction]
            },
            remainder: true
        )
        .map { _id, meta, paired, mate1, mate2 ->
            return [meta, paired, mate1 ?: [], mate2 ?: []]
        }

    DETECT(ch_combined_junctions, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(DETECT.out.versions)

    UNIFY(DETECT.out.reads
        .join(DETECT.out.coordinates)
        .join(DETECT.out.counts)
    )
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed      = UNIFY.out.bed.map{ meta, bed -> [meta + [tool: "circtools"], bed] }
    versions = ch_versions
}
