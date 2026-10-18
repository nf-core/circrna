include { CIRI2 as DETECT } from '../../../../modules/local/ciri2'
include { GAWK  as UNIFY  } from '../../../../modules/nf-core/gawk'

workflow CIRI2 {
    take:
    ch_sam
    ch_fasta
    ch_gtf

    main:
    ch_versions = channel.empty()

    DETECT(ch_sam, ch_fasta.first(), ch_gtf.first())
    ch_versions = ch_versions.mix(DETECT.out.versions)

    UNIFY(
        DETECT.out.txt.map { meta, txt ->
            [meta + [tool: "ciri2"], txt, "ciri2.bed"]
        },
        [],
        false,
    )

    emit:
    bed       = UNIFY.out.output
    ciri2_txt = DETECT.out.txt
    ciri2_sam = ch_sam
    versions  = ch_versions
}
