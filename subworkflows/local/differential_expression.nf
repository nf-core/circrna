include { HISAT2_ALIGN              } from '../../modules/nf-core/hisat2/align/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { STRINGTIE_STRINGTIE       } from '../../modules/nf-core/stringtie/stringtie/main'

workflow DIFFERENTIAL_EXPRESSION {

    take:
    reads
    gtf
    fasta
    hisat2_index

    main:
    ch_versions = Channel.empty()

    //
    // LINEAR RNA ALIGNEMT WORKFLOW:
    //

    HISAT2_EXTRACTSPLICESITES( gtf )
    HISAT2_ALIGN( reads, gtf, HISAT2_EXTRACTSPLICESITES.out.txt )
    STRINGTIE( HISAT2_ALIGN.out.bam, gtf )

    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(STRINGTIE.out.versions)

    emit:
    foo = "foo"
    versions = ch_versions
}
