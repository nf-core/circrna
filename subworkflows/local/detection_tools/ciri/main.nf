include { CIRIQUANT                     } from '../../../../modules/local/ciriquant/ciriquant'
include { CIRI_CIRI2 as CIRI2           } from '../../../../modules/local/ciri/ciri2'
include { GAWK as UNIFY                 } from '../../../../modules/nf-core/gawk'

workflow CIRI {
    take:
    ch_sam
    ch_fasta
    ch_gtf

    main:
    ch_versions = channel.empty()

    CIRI2(ch_sam, ch_fasta, ch_gtf)
    ch_versions = ch_versions.mix(CIRI2.out.versions)

    UNIFY(
        CIRI2.out.txt.map { meta, txt ->
            [meta + [tool: "ciri"], txt, "ciri.bed"]
        },
        [],
        false,
    )

    emit:
    bed                = UNIFY.out.output
    ciri_txt           = CIRI2.out.txt
    ciri_sam           = ch_sam
    versions           = ch_versions
}
