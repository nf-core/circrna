include { CIRCRNA_FINDER as MAIN } from '../../../modules/local/circrna_finder'
include { GAWK as UNIFY          } from '../../../modules/nf-core/gawk'

workflow CIRCRNA_FINDER {
    take:
    fasta
    star_sam
    star_junctions
    star_tab
    
    main:
    ch_versions = Channel.empty()

    ch_joined = star_sam.join(star_junctions).join(star_tab)
        .map{ meta, sam, junction, tab -> 
            [ meta + [tool: "circrna_finder"], [sam, junction, tab] ] }
    
    MAIN( ch_joined )
    UNIFY( MAIN.out.results, [] )

    ch_versions = ch_versions.mix(MAIN.out.versions)
    
    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}