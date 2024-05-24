include { CIRCEXPLORER2_REFERENCE as REFERENCE } from '../../../modules/local/circexplorer2/reference'
include { MAPSPLICE_ALIGN as ALIGN             } from '../../../modules/local/mapsplice/align'
include { CIRCEXPLORER2_PARSE as PARSE         } from '../../../modules/nf-core/circexplorer2/parse'
include { CIRCEXPLORER2_ANNOTATE as ANNOTATE   } from '../../../modules/nf-core/circexplorer2/annotate'
include { GAWK as UNIFY                        } from '../../../modules/nf-core/gawk'

workflow MAPSPLICE {
    take:
    reads
    gtf
    fasta
    bowtie_index
    chromosomes
    star_junctions
    
    main:
    ch_versions = Channel.empty()

    REFERENCE( gtf )
    ALIGN( reads, bowtie_index, chromosomes, gtf )
    PARSE( ALIGN.out.raw_fusions )
    ANNOTATE( PARSE.out.junction, fasta, REFERENCE.out.txt )
    UNIFY( ANNOTATE.out.txt.map{ meta, txt ->
        [ meta + [tool: "mapsplice"], txt ] }, [] )

    ch_versions = ch_versions.mix(REFERENCE.out.versions)
    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)
    
    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}