include { CIRCEXPLORER2_REFERENCE as REFERENCE } from '../../../modules/local/circexplorer2/reference'
include { MAPSPLICE_ALIGN as ALIGN             } from '../../../modules/local/mapsplice/align'
include { CIRCEXPLORER2_PARSE as PARSE         } from '../../../modules/nf-core/circexplorer2/parse'
include { CIRCEXPLORER2_ANNOTATE as ANNOTATE   } from '../../../modules/nf-core/circexplorer2/annotate'
include { CIRCEXPLORER2_FILTER as FILTER       } from '../../../modules/local/circexplorer2/filter'


workflow MAPSPLICE {
    take:
    reads
    gtf
    fasta
    bowtie_index
    chromosomes
    star_junctions
    bsj_reads
    
    main:
    ch_versions = Channel.empty()

    REFERENCE( gtf )
    ALIGN( reads, bowtie_index, chromosomes, gtf )
    PARSE( ALIGN.out.raw_fusions )
    ANNOTATE( PARSE.out.junction, fasta, REFERENCE.out.txt )
    FILTER( ANNOTATE.out.txt.map{ meta, txt ->
        [ meta + [tool: "mapsplice"], txt ] }, bsj_reads )

    ch_versions = ch_versions.mix(REFERENCE.out.versions)
    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(FILTER.out.versions)
    
    emit:
    matrix = FILTER.out.matrix
    results = FILTER.out.results

    versions = ch_versions
}