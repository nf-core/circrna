include { CIRCEXPLORER2_PARSE as PARSE         } from '../../../modules/nf-core/circexplorer2/parse'
include { GAWK as UNIFY                        } from '../../../modules/nf-core/gawk'

workflow CIRCEXPLORER2 {
    take:
    fasta
    star_junctions
    circexplorer2_index

    main:
    ch_versions = Channel.empty()

    PARSE( star_junctions )
    ch_versions = ch_versions.mix(PARSE.out.versions)

    UNIFY( PARSE.out.junction
        .map{ meta, txt -> [ meta + [tool: "circexplorer2"], txt ] }, [], false )
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed   = UNIFY.out.output

    versions = ch_versions
}
