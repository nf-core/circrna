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
    circexplorer2_index

    main:
    ch_versions = Channel.empty()

    ALIGN( reads, bowtie_index, chromosomes, gtf )
    PARSE( ALIGN.out.raw_fusions )
    ANNOTATE( PARSE.out.junction, fasta, circexplorer2_index )
    UNIFY( ANNOTATE.out.txt.map{ meta, txt ->
        [ meta + [tool: "mapsplice"], txt ] }, [] )

    ch_versions = ch_versions.mix(ALIGN.out.versions)
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)
    ch_versions = ch_versions.mix(UNIFY.out.versions)

    emit:
    bed = UNIFY.out.output

    versions = ch_versions
}
