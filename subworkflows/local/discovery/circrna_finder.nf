include { CIRCRNA_FINDER_FILTER as FILTER } from '../../../modules/local/circrna_finder/filter'

workflow CIRCRNA_FINDER {
    take:
    fasta
    star_sam
    star_junctions
    star_tab
    bsj_reads
    
    main:
    ch_versions = Channel.empty()

    ch_joined = star_sam.join(star_junctions).join(star_tab)
        .map{ meta, sam, junction, tab -> 
            [ meta + [tool: "circrna_finder"], sam, junction, tab ] }
    
    FILTER( ch_joined, fasta, bsj_reads )

    ch_versions = ch_versions.mix(FILTER.out.versions)
    
    emit:
    matrix = FILTER.out.matrix
    results = FILTER.out.results

    versions = ch_versions
}