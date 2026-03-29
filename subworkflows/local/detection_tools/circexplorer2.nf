include { CIRCEXPLORER2_PARSE as PARSE         } from '../../../modules/nf-core/circexplorer2/parse'
include { GAWK as UNIFY                        } from '../../../modules/nf-core/gawk'

workflow CIRCEXPLORER2 {
    take:
    _fasta
    star_junctions
    _circexplorer2_index

    main:
    ch_versions = channel.empty()

    PARSE( star_junctions )

    UNIFY( PARSE.out.junction
        .map{ meta, txt -> [ meta + [tool: "circexplorer2"], txt ] }, [], false )

    emit:
    bed   = UNIFY.out.output

    versions = ch_versions
}
