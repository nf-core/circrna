include { HISAT2_ALIGN              } from '../../modules/nf-core/hisat2/align/main'
include { STRINGTIE_STRINGTIE       } from '../../modules/nf-core/stringtie/stringtie/main'

workflow DIFFERENTIAL_EXPRESSION {

    take:
    reads
    gtf
    fasta
    hisat2_index
    splice_sites

    main:
    ch_versions = Channel.empty()

    //
    // LINEAR RNA ALIGNEMT WORKFLOW:
    //

    HISAT2_ALIGN( reads, hisat2_index, splice_sites )
    STRINGTIE_STRINGTIE( HISAT2_ALIGN.out.bam, gtf )

    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

    emit:
    versions = ch_versions
}
