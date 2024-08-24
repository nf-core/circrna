include { GNU_SORT as COMBINE_TRANSCRIPTOME_GTFS } from '../../modules/nf-core/gnu/sort'
include { GAWK as EXCLUDE_OVERLONG_TRANSCRIPTS   } from '../../modules/nf-core/gawk'
include { GFFREAD as TRANSCRIPTOME               } from '../../modules/local/gffread'

workflow COMBINE_TRANSCRIPTOMES {
    take:
    ch_genome_fasta
    ch_genome_gtf
    ch_circ_gtf
    
    main:
    ch_versions = Channel.empty()
    
    COMBINE_TRANSCRIPTOME_GTFS(
        ch_genome_gtf.mix(ch_circ_gtf).map{meta, gtf -> gtf}.collect().map{[[id: "transcriptome"], it]},
    )
    ch_versions = ch_versions.mix(COMBINE_TRANSCRIPTOME_GTFS.out.versions)

    EXCLUDE_OVERLONG_TRANSCRIPTS(
        COMBINE_TRANSCRIPTOME_GTFS.out.sorted, []
    )
    ch_versions = ch_versions.mix(EXCLUDE_OVERLONG_TRANSCRIPTS.out.versions)

    TRANSCRIPTOME(EXCLUDE_OVERLONG_TRANSCRIPTS.out.output, ch_genome_fasta)
    ch_versions = ch_versions.mix(TRANSCRIPTOME.out.versions)


    emit:
    fasta = TRANSCRIPTOME.out.output
    gtf   = EXCLUDE_OVERLONG_TRANSCRIPTS.out.output

    versions = ch_versions
}