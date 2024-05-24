include { CIRCEXPLORER2_REFERENCE as REFERENCE } from '../../../modules/local/circexplorer2/reference'
include { CIRCEXPLORER2_PARSE as PARSE         } from '../../../modules/nf-core/circexplorer2/parse'
include { CIRCEXPLORER2_ANNOTATE as ANNOTATE   } from '../../../modules/nf-core/circexplorer2/annotate'
include { GAWK as UNIFY                        } from '../../../modules/nf-core/gawk'

workflow CIRCEXPLORER2 {
    take:
    gtf
    fasta
    star_junctions
    
    main:
    ch_versions = Channel.empty()
    
    REFERENCE( gtf )
    PARSE( star_junctions )
    ANNOTATE( PARSE.out.junction, fasta, REFERENCE.out.txt )
    UNIFY( ANNOTATE.out.txt
        .map{ meta, txt -> [ meta + [tool: "circexplorer2"], txt ] }, [] )

    ch_versions = ch_versions.mix(REFERENCE.out.versions)
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed   = UNIFY.out.output

    versions = ch_versions
}