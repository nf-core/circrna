include { MAPSPLICE_ALIGN as ALIGN             } from '../../../modules/local/mapsplice/align'
include { CIRCEXPLORER2_PARSE as PARSE         } from '../../../modules/nf-core/circexplorer2/parse'
include { GAWK as UNIFY                        } from '../../../modules/nf-core/gawk'

workflow MAPSPLICE {
    take:
    reads
    gtf
    fasta
    bowtie_index
    chromosomes
    star_junctions
    circexplorer2_index

    main:
    ch_versions = Channel.empty()

    ALIGN( reads, bowtie_index.map{ _meta, index -> index}, chromosomes, gtf )
    PARSE( ALIGN.out.raw_fusions )
    UNIFY( PARSE.out.junction.map{ meta, bed ->
        [ meta + [tool: "mapsplice"], bed ] }, [], false )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}
